import sys
import os
from argparse import ArgumentParser
from sklearn.model_selection import train_test_split

ALPHABET = ["A", "C", "G", "T"]
CLS_TOKEN = "[CLS]"
SEP_TOKEN = "[SEP]"
UNK_TOKEN = "[UNK]"


def preprocess_gencode(file_name):
    with open(file_name) as f:
        data = f.readlines()
    chr_list = []
    cur_seq = ""
    for item in data:
        if item[0:4] == ">chr":
            print("Reading {}...".format(item.strip()))
            if len(cur_seq) > 0:
                chr_list.append(cur_seq)
            cur_seq = ""
        else:
            cur_seq += item.strip().upper()
    return chr_list


def create_data_split(chr_list, max_len):
    agg_list = []
    for chr in chr_list:
        agg_list += [chr[i:i+max_len] for i in range(0, len(chr), max_len)]
    train_seq, test_seq = train_test_split(agg_list, test_size=0.2, random_state=1)
    train_seq, val_seq = train_test_split(train_seq, test_size=0.25, random_state=1)
    return {"train": train_seq, "val": val_seq, "test": test_seq}


def seq2kmer(seq, k):
    kmer = [seq[i:i+k] for i in range(len(seq)+1-k)]
    return kmer


def write_to_file(split_dict, output_dir, k):
    for key in split_dict:
        cur_bucket = split_dict[key]
        with open(os.path.join(output_dir, "{}_seq.txt".format(key)), "w") as f:
            for seq in cur_bucket:
                f.write(CLS_TOKEN + " ")
                kmer_list = seq2kmer(seq, k)
                for kmer in kmer_list:
                    if all(nucl in ALPHABET for nucl in kmer):
                        f.write(kmer + " ")
                    else:
                        f.write(UNK_TOKEN + " ")
                f.write(SEP_TOKEN + "\n")


if __name__ == "__main__":
    parser = ArgumentParser(description="Preprocessing script for GENCODE GRCh38.p13 data.")
    parser.add_argument("-i", "--input_file", help="Input genome file.")
    parser.add_argument("-o", "--output_dir", help="Output directory.")
    parser.add_argument("-l", "--max_len", help="Maximum sequence length.", type=int, default=1024)
    parser.add_argument("-k", "--k_mer", help="k value for k-mer.", type=int, default=3)
    args = parser.parse_args(sys.argv[1:])

    chr_list = preprocess_gencode(args.input_file)
    split_dict = create_data_split(chr_list, max_len=args.max_len)
    write_to_file(split_dict, args.output_dir, args.k_mer)
