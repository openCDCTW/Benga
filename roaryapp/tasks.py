from roaryapp.celeryapp import app
import os
import os.path


def format_cmd(prog, flags, arg):
    a = list(map(lambda x: x[0] + "=" + str(x[1]) if x[1] != "" else x[0], flags))
    flag = " ".join(a)
    return " ".join([prog, flag, arg])


@app.task
def async_roary(uploads, ident_min, threads):
    # write uploaded file
    for filename, file in uploads:
        with open(os.path.join("/input", filename), "w") as f:
            f.write(file)

    # format command and execute
    output_dir = "/output"
    args = list()
    args.append(("-p", threads))
    args.append(("-i", ident_min))
    args.append(("-f", os.path.join(output_dir, "roaryapp")))
    cmd = format_cmd("roaryapp", args, os.path.join("/input", "*.gff"))
    os.system(cmd)

    # read results for return
    target_file = "gene_presence_absence.csv"
    with open(os.path.join(output_dir, "roaryapp", target_file), "r") as f:
        file = f.read()

    # clean up folder and return
    os.system("rm -r {}".format(output_dir))
    return target_file, file


