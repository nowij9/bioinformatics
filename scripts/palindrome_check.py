def palindrome_check(sequencefile):

    with open(sequencefile, "r") as fi:
        lines = fi.readlines()
        seq = ""
        for line in lines:
            if not line.startswith(">"):
                seq += line

            seq = seq.lower()
            seq = seq.replace("\n", "")
            seq = seq.replace(" ", "")
            print(seq)

        S = len(seq)

        for i in range(S // 2):
            if seq[i] != seq[S - i - 1]:
                print("Input sequence is not palindromic.")
                return

        print("Input sequence is palindromic.")


sequencefile = input()
palindrome_check(sequencefile)
