import sys
import os
import subprocess
import queue
import re
import termios
import tty
import time
from unittest import TestLoader
from pathlib import Path
from multiprocessing import cpu_count
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

OK_RE   = re.compile(r'(ok|expected failure)$')
FAIL_RE = re.compile(r'FAIL$')

# ANSI color codes
GREEN = '\033[32m'
BLUE = "\033[34m"
RED = "\033[31m"
RESET = '\033[0m'

def supports_color():
    """Returns True if the running terminal supports colors."""
    return sys.stdout.isatty() and os.getenv('TERM', '') != 'dumb'


def run_test_file(path, render=False):
    '''Run a single test file using unittest'''
    command = [sys.executable, '-m', 'unittest', '-v', str(path)]
    env = os.environ.copy()
    if render:
        env['RENDER_RMD'] = 'TRUE'
        env['RENDER_QMD'] = 'TRUE'
    try:
        result = subprocess.run(
            command, text=True, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        # count sucessful tests
        passed = sum(1 for line in result.stdout.splitlines() if OK_RE.search(line))
        failed = sum(1 for line in result.stdout.splitlines() if FAIL_RE.search(line))
        return (path.name, passed, failed, result.returncode, result.stdout)
    except Exception as e:
        return (path.name, 0, 0, 1, f'ERROR: {e}')


def get_plural(n, singular, plural):
    return f'{n} {singular if n == 1 else plural}'


def print_final_summary(n_files, n_files_passed, n_files_failed,
                        n_tests, n_tests_passed, n_tests_failed,
                        start_time):
    """Print the final summary of test results."""
    end_time = time.time()
    elapsed_time = end_time - start_time
    if elapsed_time < 60:
        time_str = f'{elapsed_time:.2f} seconds'
    else:
        mins = int(elapsed_time // 60)
        secs = int(elapsed_time % 60)
        msecs = int((elapsed_time - mins * 60 - secs) * 1000)
        time_str = f'{mins}:{secs:02d}.{msecs:02d}'

    color = supports_color()

    print(f'\nSummary\n{"-" * 80}')
    print(f'Ran {n_tests} tests in {time_str}')
    if n_files_failed > 0:
        print(f"{RED if color else ''}{get_plural(n_tests_failed, 'test', 'tests')} failed in {get_plural(n_files_failed, 'file', 'files')}.{RESET if color else ''}")
    print(f"{GREEN if color else ''}{n_files_passed} of {get_plural(n_files, 'test file', 'test files')} passed.{RESET if color else ''}")
    print(f"{GREEN if color else ''}{n_tests_passed} of {get_plural(n_tests, 'test', 'tests')} passed.{RESET if color else ''}")


def update_line(display_row, text):
    """Move the cursor to `display_row`, clear the line, write text, then restore."""
    # ANSI: save cursor, move absolute, clear line, write, restore cursor
    sys.stdout.write(f"\0337\033[{display_row};0H\033[2K{text}\0338")
    sys.stdout.flush()


def run_test_file_stream(path, q, display_row, render=False):
    """Run one test file, stream its output and post progress to a queue."""
    cmd = [sys.executable, '-m', 'unittest', '-v', str(path)]
    env = os.environ.copy()
    if render:
        env['RENDER_RMD'] = 'TRUE'
        env['RENDER_QMD'] = 'TRUE'

    p = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, env=env
    )

    passed = failed = 0
    for line in p.stdout:
        if OK_RE.search(line):
            passed += 1
        elif FAIL_RE.search(line):
            failed += 1
        q.put((path.name, display_row, passed, failed, False))

    # final “done” message
    code = p.wait()
    q.put((path.name, display_row, passed, failed, True, code))


def _print_test_file_summary(file_i, n_files, name, passed, failed, n_tests, verbose=False):
    if verbose:
        print('-' * 80)
    icon = '✅' if failed == 0 else '❌'
    len_n_files = len(str(n_files))
    test_str = f'({passed}/{n_tests})'
    print(f'{icon} File {str(file_i).ljust(len_n_files)} of {n_files} {test_str.ljust(9)} {name}')
    if verbose:
        print('-' * 80)
    sys.stdout.flush()


def _run_non_interactive(test_files, file_test_counts, n_cores, verbose=False, **kwargs):
    n_files = len(test_files)
    n_tests = sum(file_test_counts.values())

    print(f'Running {get_plural(n_files, 'test file', 'test files')} on {get_plural(n_cores, 'core', 'cores')}...')

    n_files_passed = n_files_failed = 0
    n_tests_passed = n_tests_failed = 0
    start_time = time.time()
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        futures = {
            executor.submit(run_test_file, f, **kwargs): f
            for f in test_files
        }
        for i, future in enumerate(as_completed(futures)):
            name, n_passed, n_failed, code, output = future.result()
            _print_test_file_summary(i + 1, n_files, name, n_passed, n_failed,
                                     file_test_counts[name], verbose=verbose)
            if code == 0:
                if verbose:
                    print(output)
                n_files_passed += 1
                n_tests_passed += n_passed
            else:
                print(output)
                n_files_failed += 1
                n_tests_passed += n_passed
                n_tests_failed += n_failed

    print_final_summary(
        n_files, n_files_passed, n_files_failed,
        n_tests, n_tests_passed, n_tests_failed,
        start_time
    )

    return 1 if n_files_failed else 0


def _get_cursor_position():
    """Ask the terminal where the cursor currently is (row, col)."""
    fd = sys.stdin.fileno()
    old = termios.tcgetattr(fd)
    try:
        tty.setraw(fd)
        sys.stdout.write('\033[6n')
        sys.stdout.flush()
        resp = ''
        while True:
            ch = sys.stdin.read(1)
            resp += ch
            if ch == 'R':
                break
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old)
    m = re.match(r'\x1b\[(\d+);(\d+)R', resp)
    if m:
        return int(m.group(1)), int(m.group(2))
    return 1, 1    # fall back if we can’t parse


def get_row_text(test_name, rjust, n_tests, started=True, n_passed=0, n_failed=0, finished=False):
    color_support = supports_color()
    ret = f'{test_name.rjust(rjust)}'

    if started:
        ret += ' ['
        if color_support:
            ret += GREEN if finished and n_failed == 0 else RED if finished else BLUE
        ret += f'{str(round(n_passed/n_tests * 100)).rjust(3)}%'
        ret += RESET if color_support else ''
        ret += f'] {n_passed} of {n_tests}'
    else:
        ret += ' -'

    if finished:
        if n_failed > 0:
            ret += ' ❌' if color_support else 'X'
        else:
            ret += GREEN if color_support else ''
            ret += ' ✔'
            ret += RESET if color_support else ''

    return ret


def _run_interactive(test_files, file_test_counts, n_cores, **kwargs):

    # 1) reserve enough lines so nothing scrolls off
    n_test_files = len(test_files)
    lines_to_reserve = n_test_files + 2
    for _ in range(lines_to_reserve):
        print()

    # 2) determine cursor position and compute where to print the header
    cur_row, _ = _get_cursor_position()

    # 3) print the header
    header_row = cur_row - (lines_to_reserve - 1)
    longest_name_len = max(20, max(len(f.name) for f in test_files))
    sys.stdout.write(f"\033[{header_row};0H")   # move to header row
    print(f'Running {get_plural(n_test_files, 'file', 'files')} on {get_plural(n_cores, 'core', 'cores')}...')

    # 4) print initial per-file lines and remember their absolute rows
    rows = {}
    for idx, f in enumerate(test_files, start=1):
        display_row = header_row + idx
        rows[f.name] = display_row
        sys.stdout.write(f"\033[{display_row};0H")   # move to the file’s row
        print(get_row_text(os.path.splitext(f.name)[0], longest_name_len, file_test_counts[f.name], started=False))
    sys.stdout.flush()

    q = queue.Queue()
    finished = n_files_passed = n_files_failed = 0
    n_tests_passed = n_tests_failed = 0

    # 5) launch threads to stream each test file
    start_time = time.time()
    with ThreadPoolExecutor(max_workers=n_cores) as executor:
        for f in test_files:
            executor.submit(
                run_test_file_stream,
                f, q, rows[f.name], kwargs.get('render', False)
            )

        # 6) process the queue until all are done
        while finished < n_test_files:
            name, row, p_cnt, f_cnt, done, *rest = q.get()
            update_line(row, get_row_text(
                os.path.splitext(name)[0], longest_name_len, file_test_counts[name],
                n_passed=p_cnt, n_failed=f_cnt, finished=done
            ))
            if done:
                code = rest[0]
                if code == 0:
                    n_files_passed += 1
                    n_tests_passed += p_cnt
                else:
                    n_files_failed += 1
                    n_tests_passed += p_cnt
                    n_tests_failed += f_cnt
                finished += 1

    # 7) final summary below everything
    end_row = header_row + n_test_files
    sys.stdout.write(f"\033[{end_row};0H")
    print('\n')
    print_final_summary(
        n_test_files, n_files_passed, n_files_failed,
        sum(file_test_counts.values()), n_tests_passed, n_tests_failed,
        start_time
    )
    sys.exit(1 if n_files_failed else 0)


def main(test_paths, max_workers=None, verbose=False, **kwargs):
    # Determine which test files to run
    if not test_paths:
        test_files = list(Path.cwd().glob('test_*.py'))
    else:
        test_files = [Path(p) for p in test_paths]
    if not test_files:
        print('No test files found to run.')
        sys.exit(1)

    n_cores = min(cpu_count() if max_workers is None else max_workers, len(test_files))

    # Count invidual tests in each file
    loader = TestLoader()
    file_test_counts = {}
    for f in test_files:
        suite = loader.discover(start_dir=str(f.parent), pattern=f.name)
        file_test_counts[f.name] = suite.countTestCases()

    # Fallback to original behavior if verbose or non-interactive
    if verbose or not sys.stdout.isatty():
        return _run_non_interactive(test_files, file_test_counts, n_cores, verbose=verbose, **kwargs)

    return _run_interactive(test_files, file_test_counts, n_cores, **kwargs)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Run unittest files in parallel.')
    parser.add_argument(
        '--render', '-r', action='store_true', default=False,
        help='Pass RENDER_RMD and RENDER_QMD to tests'
    )
    parser.add_argument(
        '--jobs', '-j', type=int,
        help='Number of parallel processes (default: CPU cores)'
    )
    parser.add_argument(
        '--verbose', '-v', action='store_true', default=False,
        help='Show full output for each file'
    )
    parser.add_argument(
        'tests', nargs='*',
        help='test_*.py files to run. By default, all test files in the current directory are run.'
    )
    args = parser.parse_args()
    main(
        test_paths=args.tests,
        max_workers=args.jobs,
        render=args.render,
        verbose=args.verbose
    )