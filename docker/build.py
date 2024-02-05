
from subprocess import run as run_subprocess

def get_git_version():
    #check whether git is installed
    result = run_subprocess(['which','git'], capture_output=True)
    if result.returncode != 0:
        return'GIT_NOT_FOUND'

    result = run_subprocess(['git', 'status', '--porcelain'], capture_output=True)
    is_dirty  = result.stdout.decode('utf-8').strip() != ''

    result = run_subprocess(['git', 'log', '-1', '--format=%cd', '--date=local'], capture_output=True)
    git_date = result.stdout.decode('utf-8').strip()

    result = run_subprocess(['git', 'rev-parse', '--verify', 'HEAD'], capture_output=True)
    git_hash = result.stdout.decode('utf-8').strip()

    return f'git_hash:"{git_hash}", last_commit:"{git_date}", uncommited_changes:{is_dirty}'

if __name__ == '__main__':
    print(get_git_version())

