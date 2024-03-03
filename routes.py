from flask import (
    Flask,
    render_template,
    redirect,
    flash,
    url_for,
    session,
    request
)

from datetime import timedelta
from sqlalchemy.exc import (
    IntegrityError,
    DataError,
    DatabaseError,
    InterfaceError,
    InvalidRequestError,
)
from werkzeug.routing import BuildError

from flask_bcrypt import Bcrypt, generate_password_hash, check_password_hash

from flask_login import (
    UserMixin,
    login_user,
    LoginManager,
    current_user,
    logout_user,
    login_required,
)

import models
from app import create_app, db, login_manager, bcrypt
from models import User
from forms import login_form, register_form
from models import fetch_record

app = create_app()


# Keeps track of current user object
@login_manager.user_loader
def load_user(user_id):
    return User.query.get(int(user_id))


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


# Clears session storage for debugging DELETE IN FINAL VERSION!!!!!!
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
    users = User.query.all()
    return render_template('admin.html', title="Admin", users=users)


# # Register
# @app.route('/register.html/', methods=("GET", "POST"), strict_slashes=False)
# def register():
#     return render_template('register.html')


# # Login
# @app.route('/login.html/', methods=("GET", "POST"), strict_slashes=False)
# def login():
#     return render_template("login.html")

# Login route
@app.route("/login/", methods=("GET", "POST"), strict_slashes=False)
def login():
    form = login_form()

    if form.validate_on_submit():
        try:
            user = User.query.filter_by(email=form.email.data).first()
            if check_password_hash(user.password, form.password.data):
                login_user(user)
                return redirect(url_for('index'))
            else:
                flash("Invalid username or password!", "danger")
        except Exception as e:
            flash(e, "danger")

    return render_template("auth.html", form=form)


# Register route
@app.route("/register/", methods=("GET", "POST"), strict_slashes=False)
def register():
    form = register_form()
    if form.validate_on_submit():
        try:
            email = form.email.data
            password = form.password.data
            username = form.username.data

            new_user = User(
                username=username,
                email=email,
                password=bcrypt.generate_password_hash(password),
            )

            db.session.add(new_user)
            db.session.commit()
            flash(f"Account created", "success")
            return redirect(url_for("login"))

        except InvalidRequestError:
            db.session.rollback()
            flash(f"Something went wrong!", "danger")

    return render_template("auth.html", form=form)


@app.route("/logout")
@login_required
def logout():
    logout_user()
    return redirect(url_for('login'))


# Routes for the db ----------------------------------------------------------------------------------------------------

@app.route("/delete_user/<int:user_id>", methods=["DELETE"])
def deleteUser(user_id):
    try:
        user = User.query.get(user_id)
        if user:
            db.session.delete(user)
            db.session.commit()
            return f"User ID {user_id} deleted"
        else:
            return f"User ID {user_id} not found", 404
    except Exception as e:
        return f"Error deleting user: {str(e)}", 500


@app.route("/change_username/<int:user_id>", methods=["PUT"])
def changeUsername(user_id):
    try:
        data = request.get_json()  # Retrieve data from request body
        new_username = data.get('new_username')  # Get the new username

        user = User.query.get(user_id)
        if user:
            user.username = new_username  # Update the username
            db.session.commit()
            return f"User ID {user_id} username changed to {new_username} successfully"
        else:
            return f"User ID {user_id} not found", 404
    except Exception as e:
        return f"Error changing username: {str(e)}", 500


# !!!!!!!!!!!!!!!!!!!!!!!!!!!! MAKE SURE TO HASH PASSWORD w BCRYPT!!!!!!!!!!!!!!!!!
@app.route("/change_password/<int:user_id>", methods=["PUT"])
def changePassword(user_id):
    try:
        data = request.get_json()  # Retrieve data from request body
        new_password = data.get('new_password')

        user = User.query.get(user_id)
        if user:
            user.password = new_password
            db.session.commit()
            return f"User ID {user_id} password changed successfully"
        else:
            return f"User ID {user_id} not found", 404
    except Exception as e:
        return f"Error changing password: {str(e)}", 500


if __name__ == "__main__":
    app.run(debug=True)
