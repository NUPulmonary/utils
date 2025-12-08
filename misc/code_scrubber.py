import os
import re
import shutil
import csv

# ---------------------------------------------
# CONFIG
# ---------------------------------------------
VALID_EXTENSIONS = {".rmd", ".qmd", ".Rmd", ".Qmd", ".sh", ".bash", ".zsh", ".r", ".R",".py", ".ipynb", ".txt", ".md"}
MARKDOWN_EXTENSIONS = {".md", ".rmd", ".Rmd", ".qmd", ".Qmd"}

# Output directory for sanitized files
OUTPUT_DIR = "cleaned_output"

# Log file
LOG_FILE = "redaction_log.tsv"

# -------- Regex patterns --------
PATTERNS = {
    "FILE_PATH": re.compile(r"[\"\']*(/projects/[bp]\d{4,5}/|~/|/mnt/|/var/|/home/[a-zA-z]{3}\d{4})\S*[\"\']*"),
    "DATE_YYMMDD": re.compile(r"\b[0-9]{6}_"),
    #autopsy ID, MRNs, SSNs, phone numbers, general long numeric IDs
    "IDENTIFIER": re.compile(r"\b(A|NMA)\d{2}-?\d+.*\b|\d{3}-\d{3}-\d{4}|\d{3}-\d{2}-\d{4}|\d{7,128}"),
    "EMAIL_URL": re.compile(r"\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.(?:com|org|edu|gov|net|mil|io|co)\b"),
}


# ---------------------------------------------
# Extract only code chunks
# ---------------------------------------------
def extract_code_chunks(text):
    """
    Retains ONLY fenced code blocks.
    """
    chunks = []
    in_chunk = False
    buffer = []

    for line in text.splitlines():
        if re.match(r"^```", line):
            if in_chunk:
                buffer.append(line)
                chunks.append("\n".join(buffer))
                buffer = []
                in_chunk = False
            else:
                buffer.append(line)
                in_chunk = True
        else:
            if in_chunk:
                buffer.append(line)

    return "\n\n".join(chunks)


# ---------------------------------------------
# Sanitize text and log redactions
# ---------------------------------------------
def sanitize_text(text, rel_path, log_rows):
    """
    Applies all patterns and logs every redaction.
    """
    for label, pattern in PATTERNS.items():
        for match in pattern.finditer(text):
            snippet = match.group(0)
            # Truncate logged text if too long
            logged_snippet = snippet if len(snippet) <= 120 else snippet[:117] + "..."

            log_rows.append({
                "file": rel_path,
                "pattern": label,
                "matched_text": logged_snippet
            })

        text = pattern.sub(f"[REDACTED_{label}]", text)

    return text


# ---------------------------------------------
# Remove leading YYMMDD_ from filenames
# ---------------------------------------------
def sanitize_filename(fn):
    return re.sub(r"^\d{6}_", "", fn)


# ---------------------------------------------
# Process a single file
# ---------------------------------------------
def process_file(src_path, dst_path, rel_path, log_rows):
    ext = os.path.splitext(src_path)[1]
    if ext not in VALID_EXTENSIONS:
        return

    with open(src_path, "r", encoding="utf-8") as f:
        original = f.read()

    # Keep only code chunks
    code_only = original
    if(ext in MARKDOWN_EXTENSIONS):
        code_only = extract_code_chunks(original)

    # Apply sanitization
    cleaned = sanitize_text(code_only, rel_path, log_rows)

    # Write sanitized output
    os.makedirs(os.path.dirname(dst_path), exist_ok=True)
    with open(dst_path, "w", encoding="utf-8") as f:
        f.write(cleaned)


# ---------------------------------------------
# Walk directory recursively
# ---------------------------------------------
def process_directory(root):
    log_rows = []

    for dirpath, _, files in os.walk(root):
        for fn in files:
            src_path = os.path.join(dirpath, fn)

            # Compute relative path for reporting
            rel_path = os.path.relpath(src_path, root)

            # Sanitize filename
            clean_fn = sanitize_filename(fn)

            dst_path = os.path.join(
                OUTPUT_DIR,
                os.path.relpath(dirpath, root),
                clean_fn
            )

            process_file(src_path, dst_path, rel_path, log_rows)

    # Write log
    with open(LOG_FILE, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["file", "pattern", "matched_text"], delimiter="\t")
        writer.writeheader()
        writer.writerows(log_rows)

    print(f"Sanitization complete.\nCleaned files in: {OUTPUT_DIR}\nRedaction log: {LOG_FILE}")


# ---------------------------------------------
# RUN
# ---------------------------------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Redact and sanitize codebase for publication.")
    parser.add_argument("directory", help="Root directory to sanitize")
    args = parser.parse_args()

    process_directory(args.directory)