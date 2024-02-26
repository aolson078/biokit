from flask import (
    render_template,
    redirect,
    url_for,
)

from app import create_app

app = create_app()


# Home
@app.route('/', methods=("GET", "POST"), strict_slashes=False)
def index():
    return render_template('index.html', title="Home")

@app.route("/index.html")
def home():
    return redirect(url_for("index"))

@app.route("/index")
def home1():
    return redirect(url_for("index"))

# Manager
@app.route('/manager/', methods=("GET", "POST"), strict_slashes=False)
def manager():
    return render_template('manager.html', title="Manager")

# Admin
@app.route('/admin/', methods=("GET", "POST"), strict_slashes=False)
def admin():
    return render_template('admin.html', title="Admin")

# Register
@app.route('/register/', methods=("GET", "POST"), strict_slashes=False)
def register():
    return render_template('register.html')

# Login
@app.route('/login/', methods=("GET", "POST"), strict_slashes=False)
def login():
    return render_template("login.html")


if __name__ == "__main__":
    app.run(debug=True)