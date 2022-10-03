import argparse
import glob
import gzip
import hashlib
import os
import re

if __name__ == "__main__":
    # Setting up argparser
    parser = argparse.ArgumentParser(description="A script to create file assertions for nf-test")
    parser.add_argument(
        "test_dir",
        metavar="TEST_DIRECTORY",
        type=str,
        help="The folder containing the test outputs (usually called `.nf-test`",
    )

    args = parser.parse_args()

    test_dir = args.test_dir
    all_outputs = glob.glob(f"{test_dir}/**/output/**", recursive=True)

    tab = "\\t"

    print("assert workflow.success")

    for output in all_outputs:
        abs_path = os.path.abspath(output)
        if re.search("^.*/multiqc_data/", output) or re.search("^.*/multiqc_plots/", output):
            continue
        if os.path.isfile(abs_path):
            file_name = re.search("^.*/output/(.*)$", output).group(1)
            if re.search("^.*\.tbi$", output) or re.search("^.*\.db$", output) or re.search("^.*multiqc_report.html$", output):
                print(f'assert file("${{outputDir}}/{file_name}").exists()')
            elif re.search("^.*\.vcf.gz$", output):
                with gzip.open(abs_path, "rt") as vcf:
                    for line in vcf:
                        if re.search("^chr.*$", line):
                            print(
                                f'assert path("${{outputDir}}/{file_name}").linesGzip.contains("{tab.join(line.split()).strip()}")'
                            )
                            break
            elif re.search("^.*\.vcf$", output):
                with open(abs_path, "r") as vcf:
                    for line in vcf:
                        if re.search("^chr.*$", line):
                            print(
                                f'assert path("${{outputDir}}/{file_name}").text.contains("{tab.join(line.split()).strip()}")'
                            )
                            break
            else:
                md5sum = hashlib.md5(open(abs_path, "rb").read()).hexdigest()
                print(f'assert path("${{outputDir}}/{file_name}").md5 == "{md5sum}"')
