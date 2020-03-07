from flask import Flask
from diplom_counting import getFirstLocationAngleCoords
from jinja2 import Template

app = Flask(__name__)


@app.route("/")
def hello():
    _, res = getFirstLocationAngleCoords()

    with open('output.html', encoding="utf-8") as file_:
        template = Template(file_.read())

    return template.render(x_s=res['sat_coords']['x_s'],
                           y_s=res['sat_coords']['y_s'],
                           z_s=res['sat_coords']['z_s'],
                           t=res['time'],
                           xi=res['xi'],
                           l_x=res['angle']['l_x'],
                           l_y=res['angle']['l_y'],
                           l_z=res['angle']['l_z'])


if __name__ == "__main__":
    app.run()
