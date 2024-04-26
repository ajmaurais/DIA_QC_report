
import os
import subprocess
import inspect


def make_work_dir(work_dir, clear_dir=False):
    if not os.path.isdir(work_dir):
        if os.path.isfile(work_dir):
            raise RuntimeError('Cannot create work directory!')
        os.makedirs(work_dir)
    else:
        if clear_dir:
            for file in os.listdir(work_dir):
                os.remove(f'{work_dir}/{file}')


def run_command(command, wd, prefix=None):
    result = subprocess.run(command, cwd=wd,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            shell=False, check=False)

    prefix_path = f'{wd}/{prefix if prefix else inspect.stack()[1][3]}'
    with open(f'{prefix_path}_command.txt', 'w') as outF:
        outF.write(f"{' '.join(command)}\n")
    with open(f'{prefix_path}.stdout.txt', 'w') as outF:
        outF.write(f"{result.stdout.decode('utf-8')}\n")
    with open(f'{prefix_path}.stderr.txt', 'w') as outF:
        outF.write(f"{result.stderr.decode('utf-8')}\n")
    with open(f'{prefix_path}.rc.txt', 'w') as outF:
        outF.write(f'{str(result.returncode)}\n')
    return result


def setup_multi_db(data_dir, output_dir, group_by_gene=False):
    make_work_dir(output_dir)
    grouping = 'by_gene' if group_by_gene else 'by_protein'

    commands = [['parse_data', '--projectName=Sp3',
                 '-m', f'{data_dir}/metadata/Sp3_metadata.tsv',
                 f'{data_dir}/skyline_reports/Sp3_replicate_quality.tsv',
                 f'{data_dir}/skyline_reports/Sp3_{grouping}_precursor_quality.tsv'],
                ['parse_data', '--overwriteMode=append', '--projectName=Strap',
                 '-m', f'{data_dir}/metadata/Strap_metadata.tsv',
                 f'{data_dir}/skyline_reports/Strap_replicate_quality.tsv',
                 f'{data_dir}/skyline_reports/Strap_{grouping}_precursor_quality.tsv']]

    if os.path.isfile(f'{output_dir}/data.db3'):
        commands[0].insert(1, '--overwriteMode=overwrite')

    if group_by_gene:
        for i in range(len(commands)):
            commands[i].insert(1, '--groupBy=gene')

    results = list()
    for command in commands:
        results.append(run_command(command, output_dir))

    return results

