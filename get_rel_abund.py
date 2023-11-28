import csv
import argparse

output_columns = ["sample", "relative_abundance", "taxonomy"]
input_columns = ["sample", "filled_coverage", "taxonomy"]

def get_args():
    parser = argparse.ArgumentParser(description="Calculate relative abundance of each taxonomy for a filled SingleM profile")
    parser.add_argument("-i", "--input", help="Input file", required=True)
    parser.add_argument("-o", "--output", help="Output file", required=True)

    return parser.parse_args()

def main():
    args = get_args()
    input_file = args.input
    output_file = args.output

    with open(input_file, "r") as r:
        reader = csv.DictReader(r, delimiter="\t", fieldnames=input_columns)
        with open(output_file, "w+") as w:
            writer = csv.DictWriter(w, delimiter="\t", fieldnames=output_columns)
            writer.writeheader()
            # get first line of reader, take filled coverage as denominator
            next(reader)
            # if end of file
            if reader.line_num == sum(1 for line in open(input_file, "r")):
                print(f"Only one line in {input_file}, skipping")
                return
            first_line = next(reader)
            denominator = float(first_line["filled_coverage"])
            print(f"Total coverage {denominator} for {input_file}")
            writer.writerow({
                "sample": first_line["sample"],
                "relative_abundance": str(denominator/denominator * 100),
                "taxonomy": first_line["taxonomy"]
            })
            for line in reader:
                # take filled coverage as numerator
                numerator = float(line["filled_coverage"])
                # calculate relative abundance
                relative_abundance = numerator / denominator * 100
                # write to output file
                writer.writerow({
                    "sample": line["sample"],
                    "relative_abundance": relative_abundance,
                    "taxonomy": line["taxonomy"]
                })
    print(f"\tFinished writing to {output_file}")

if __name__ == "__main__":
    main()