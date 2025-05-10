import math

def clean_format(input_path: str, output_path: str):
    with open(input_path, "r") as f:
        lines = f.readlines()

    output_lines = []
    line_index = 0
    for line in lines:
        line = line.split('!')[0].strip()
        if not line:
            continue
        tokens = line.split()
        evaluated = []

        for idx, tok in enumerate(tokens):
            try:
                val = eval(tok, {"sqrt": math.sqrt})
                # Force integer for specific positions:
                if line_index == 0:  # First line: atom count and dimension
                    val = int(val)
                    evaluated.append(str(val))
                elif line_index >= 1 + 2:  # After lattice lines: atom block
                    if idx == 0:
                        val = int(val)  # Atomic number
                        evaluated.append(str(val))
                    else:
                        evaluated.append(f"{val:.6f}")
                else:
                    evaluated.append(f"{val:.6f}")
            except Exception:
                continue

        output_lines.append(" ".join(evaluated))
        line_index += 1

    with open(output_path, "w") as f:
        f.write("\n".join(output_lines))

    return output_path

# Run with type-controlled storage
input_file = "Input/graphene.in"
output_file = "Input/graphene_clean.in"
clean_format(input_file, output_file)
