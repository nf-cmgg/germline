import argparse
import glob
import hashlib
import os
import re

if __name__ == "__main__":
    # Setting up argparser
    parser = argparse.ArgumentParser(description="A script to create file assertions for nf-test")
    parser.add_argument('test_dir', metavar='TEST_DIRECTORY', type=str, help="The folder containing the test outputs (usually called `.nf-test`")

    args = parser.parse_args()

    test_dir = args.test_dir
    all_outputs = glob.glob(f'{test_dir}/**/output/**', recursive=True)

    for output in all_outputs:
        abs_path = os.path.abspath(output)
        if os.path.isfile(abs_path):
            md5sum = hashlib.md5(open(abs_path, 'rb').read()).hexdigest()
            file_name = re.search('^.*/output/(.*)$', output).group(1)
            print(f"assert path(\"${{outputDir}}/{file_name}\").md5 == \"{md5sum}\"")