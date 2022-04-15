
import json
import re
from IPython.display import Markdown

def make_ipynb_toc(notebook):

    toc = []

    with open(notebook, "r") as f:
            cells = json.load(f)["cells"]

    for cell in cells:
        if cell["cell_type"] == "markdown":
            for line in cell["source"]:
                match = re.search("^#+ \w+", line)
                if match:
                    level = len(line) - len(line.lstrip("#"))
                    link = line.strip(" #\n").replace(" ", "-")
                    toc.append(
                        2 * (level - 1) * " "
                        + "- ["
                        + line.strip(" #\n")
                        + "](#"
                        + link
                        + ")"
                    )
    return(Markdown('\n'.join(toc)))