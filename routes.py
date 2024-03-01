from flask import (
    render_template,
    redirect,
    url_for,
    request,
    session
)

from app import create_app
from models import fetch_record

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

@app.route('/clear_session', methods=['GET'])
def clear_session():
    session.clear()
    return redirect(url_for("employee"))



@app.route("/search", methods=["POST"], strict_slashes=False)
def search():
    # Get data from employee form
    query = request.form.get("search")
    results = fetch_record(query)
    return render_template('results.html', query=query, results=results)

# Employee
@app.route('/employee.html/', methods=("GET", "POST"), strict_slashes=False)
@app.route('/employee.html/<selected_result>', methods=("GET", "POST"), strict_slashes=False)
def employee(selected_result=None):
    if selected_result:
        if 'results' not in session:
            session['results'] = []
        # Store results in session
        session['results'].append(selected_result)
        session.modified = True
    return render_template('employee.html', title="Employee")


# Manager
@app.route('/manager.html/', methods=("GET", "POST"), strict_slashes=False)
def manager():
    return render_template('manager.html', title="Manager")


# Admin
@app.route('/admin.html/', methods=("GET", "POST"), strict_slashes=False)
def admin():
    return render_template('admin.html', title="Admin")


# Register
@app.route('/register.html/', methods=("GET", "POST"), strict_slashes=False)
def register():
    return render_template('register.html')


# Login
@app.route('/login.html/', methods=("GET", "POST"), strict_slashes=False)
def login():
    return render_template("login.html")


if __name__ == "__main__":
    app.secret_key = 'your-secret-key'
    app.run(debug=True)
