from wtforms import (
	StringField,
	PasswordField,
)

from flask_wtf import FlaskForm
from wtforms.fields.core import RadioField
from wtforms.validators import InputRequired, Length, EqualTo, Email
from wtforms import ValidationError
from models import User


class register_form(FlaskForm):
	username = StringField(
		validators=[
			InputRequired(),
			Length(1, 50, message="Please provide a valid name"),
		]
	)
	email = StringField(validators=[InputRequired(), Email(), Length(5, 100)])
	password = PasswordField(validators=[InputRequired(), Length(2, 100)])
	confirm_password = PasswordField(
		validators=[
			InputRequired(),
			Length(2, 100),
			EqualTo("password", message="Passwords must match !"),
		]
	)

	role = RadioField('Role',
	                           choices=[('employee', 'Employee - Standard user, can query and compile reports'),
	                                    ('manager', 'Manager - User with privileges to manage and print reports'),
	                                    ('admin', 'Admin - User with root privileges, can create other accounts')],

	                           validators=[
									InputRequired()
	                                        ])

	# checks if email is already registered
	def valid_email(self, email):
		if User.query.filter_by(email=email.data).first():
			raise ValidationError("Email already registered!")

	# checks if username is already registered
	def valid_username(self, username):
		if User.query.filter_by(username=username.data).first():
			raise ValidationError("Username already taken!")


class login_form(FlaskForm):
	email = StringField(validators=[InputRequired(), Email(), Length(5, 100)])
	password = PasswordField(validators=[InputRequired(), Length(2, 100)])
