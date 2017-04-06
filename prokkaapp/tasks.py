from prokkaapp.celeryapp import app
import os
import os.path


def format_cmd(prog, flags, arg):
    a = list(map(lambda x: x[0] + "=" + str(x[1]) if x[1] != "" else x[0], flags))
    flag = " ".join(a)
    return " ".join([prog, flag, arg])


@app.task
def async_prokka(filename, file):
    # write uploaded file
    with open(os.path.join("/input", filename), "w") as f:
        f.write(file)

    # format command and execute
    name, ext = filename.split(".")
    output_dir = os.path.join("/output", name)
    args = list()
    args.append(("--cpus", "2"))
    args.append(("--outdir", output_dir))
    args.append(("--prefix", name))
    cmd = format_cmd("prokkaapp", args, os.path.join("/input", filename))
    os.system(cmd)

    # read results for return
    filepath = os.path.join(output_dir, name)
    with open("{}.ffn".format(filepath), "r") as f:
        ffn_file = f.read()
    with open("{}.gff".format(filepath), "r") as f:
        gff_file = f.read()
    result = [("{}.ffn".format(name), ffn_file), ("{}.gff".format(name), gff_file)]

    # clean up folder and return
    os.system("rm -r {}".format(output_dir))
    return result
