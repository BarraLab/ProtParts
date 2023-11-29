import os

MAKEBLASTDB_EXEC = 'makeblastdb'
BLASTP_EXEC = 'blastp'
TMP_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'tmp')

if not os.path.exists(TMP_DIR):
    os.mkdir(TMP_DIR)