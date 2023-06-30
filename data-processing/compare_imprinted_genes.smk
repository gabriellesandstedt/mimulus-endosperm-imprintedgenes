rule shared_callable_genes:

til_MEGs="til_callable.txt"
caes_MEGs="caes_callable.txt"
til_MEGs1="til_callable1.txt"
caes_MEGs1="caes_callable1.txt"
shared_MEGs="shared_callable.txt"

awk '{print $1}' "$til_MEGs" > "$til_MEGs1"
awk '{print $1}' "$caes_MEGs" > "$caes_MEGs1"

grep -Fxf "$til_MEGs1" "$caes_MEGs1" > "$shared_MEGs"

# shared til MEGs + callable
# shared caes MEGs + callable

# shared til PEGs + callable
# shared caes PEGs + callable 

# shared caes and til MEGs and PEGs

# shared til megs CPM
  remove quotations from :  edgeR_result_sign_til.txt
  sed 's/"//g' edgeR_result_sign_til.txt > til_CPM.txt

  pull out MEGs CPM: 
  
# shared caes megs CPM


rule shared_MEGs:
    input:
        til_MEGs="edgeR_List_MEGs_SOP12xLVR1_LFC2.txt",
        caes_MEGs="edgeR_List_MEGs_UTC1xTWN36_LFC2.txt"
    output:
        shared_MEGs="shared_MEGs.txt"
    shell:
        """
        grep -Fxf {input.til_MEGs} {input.caes_MEGs} > {output.shared_MEGs}
        """
        
rule shared_PEGs:
    input:
        til_PEGs="edgeR_List_PEGs_SOP12xLVR1_LFC2.txt",
        caes_PEGs="edgeR_List_PEGs_UTC1xTWN36_LFC2.txt"
    output:
        shared_PEGs="shared_PEGs.txt"
    shell:
        """
        grep -Fxf {input.til_PEGs} {input.caes_PEGs} > {output.shared_PEGs}
        """

rule shared_MEGs_shared_callable_genes:

rule shared_PEGs_shared_callable_genes:
        
rule unique_til_MEGs:
    input:
        til_MEGs="file1.txt",
        caes_MEGs="file2.txt"
    output:
        unique_til_MEGs="unique_til_MEGs.txt"
    shell:
        """
        file1={input.til_MEGs}
        file2={input.caes_MEGs}

        while read -r line || [ -n "$line" ]; do
            if ! grep -q "^$line$" "$file2"; then
                echo "$line"
            fi
        done < "$file1" > {output.unique_til_MEGs}
        """

rule unique_caes_MEGs:
    input:
        caes_MEGs="file1.txt",
        til_MEGs="file2.txt"
    output:
        unique_caes_MEGs="unique_caes_MEGs.txt"
    shell:
        """
        file1={input.til_MEGs}
        file2={input.caes_MEGs}

        while read -r line || [ -n "$line" ]; do
            if ! grep -q "^$line$" "$file2"; then
                echo "$line"
            fi
        done < "$file1" > {output.unique_caes_MEGs}
        """

rule unique_til_PEGs:
    input:
        til_PEGs="file1.txt",
        caes_PEGs="file2.txt"
    output:
        unique_til_PEGs="unique_til_PEGs.txt"
    shell:
        """
        file1={input.til_PEGs}
        file2={input.caes_PEGs}

        while read -r line || [ -n "$line" ]; do
            if ! grep -q "^$line$" "$file2"; then
                echo "$line"
            fi
        done < "$file1" > {output.unique_til_PEGs}
        """

rule unique_caes_PEGs:
    input:
        caes_MEGs="file1.txt",
        til_MEGs="file2.txt"
    output:
        unique_caes_PEGs="unique_caes_PEGs.txt"
    shell:
        """
        file1={input.til_PEGs}
        file2={input.caes_PEGs}

        while read -r line || [ -n "$line" ]; do
            if ! grep -q "^$line$" "$file2"; then
                echo "$line"
            fi
        done < "$file1" > {output.unique_caes_PEGs}
        """

rule extract_normalized_counts_for_til_MEGs:
    input:
        file1 = 'file1.txt',
        file2 = 'file2.txt'
    output:
        'matching_rows.txt'
    script:
        """
        def extract_matching_rows(file1, file2):
            with open(file2, 'r') as f2:
                rows_set = set(line.strip() for line in f2)

            with open(file1, 'r') as f1:
                matching_rows = [line.strip() for line in f1 if line.strip().split()[12] in rows_set]

            return matching_rows

        matching_rows = extract_matching_rows('{input.file1}', '{input.file2}')

        with open('{output}', 'w') as outfile:
            for row in matching_rows:
                outfile.write(row + '\\n')
        """

