# BioKit Login 500 Handoff

Date: 2026-08-25

## Current outcome

BioKit is running at `http://127.0.0.1:5000`, but login with the README credentials still returns HTTP 500.

The latest confirmed server error is:

```text
sqlite3.OperationalError: no such column: user.password_hash
```

The application now reaches the shipped SQLite database and finds the `user` table, but the database schema does not match the current SQLAlchemy `User` model. The model selects `user.password_hash`; that column is absent from the shipped database.

## Confirmed failure sequence

1. The login form originally failed before authentication because WTForms could not import `email_validator`.
2. After declaring that dependency, the login route failed because `User` was not imported in the route scope.
3. After adding the import, the app used Flask's instance-relative default database and failed with `no such table: user`.
4. After pointing the default configuration at the shipped root `database.db`, the app reached the existing `user` table but failed with `no such column: user.password_hash`.

The fourth failure is the current blocker.

## Uncommitted changes

- `requirements.txt`: added `email-validator`.
- `flask_bio_app.py`: the default SQLite URI now resolves to the shipped root `database.db`.
- `flask_bio_app.py`: the login route now imports `User` locally, avoiding the circular-import issue.
- `tests/test_search_and_auth.py`: added a focused successful-login test using an isolated in-memory SQLite database.

Current working-tree state also includes an untracked `instance/` directory created when Flask opened the former instance-relative database. It has not been deleted or altered.

## Verification completed

- Docker web, worker, and Redis containers are running.
- Active web configuration was verified as `sqlite:////app/database.db`.
- `email_validator` imports successfully in the rebuilt web image.
- WTForms accepts a valid email and rejects an invalid email.
- Focused login test passed against an isolated in-memory database:

```text
pytest -q tests/test_search_and_auth.py::test_login_accepts_valid_credentials
1 passed
```

- A reversible fault injection removed the `User` import; the focused test reproduced the exact `NameError`. The import was restored and the clean test passed again.
- The live browser login has not passed. Current logs prove the shipped database/model mismatch above.

## Why this was not completed

The next change is no longer a normal source edit. It requires deciding how to reconcile an existing database schema and its stored password representation with the current model. The repository has Flask-Migrate scaffolding but no versioned migration scripts, so there is no declared migration to apply.

Changing, renaming, or backfilling authentication columns can alter persistent account data. The applicable safety rules require:

1. Explicit permission for database inspection, including the exact local database target.
2. Read-only inspection of the existing `user` schema and password-related columns.
3. A verified backup outside the migration scope and a precise rollback path.
4. A preview of the exact migration, target database, affected rows, and expected result.
5. A separate approval before executing the migration.

Those steps were not authorized, so no database mutation was attempted.

Browser automation also did not provide end-to-end proof. Computer Use could not confidently verify Brave's current URL and stopped before navigating or submitting the login form. No automated browser login was completed.

## Recommended next steps

1. Authorize a read-only schema inspection of the local file `database.db`, beginning with the equivalent of `PRAGMA table_info(user)`; do not modify data.
2. Determine the existing password column name and storage format. Do not assume it is plaintext or compatible with Werkzeug hashes.
3. Choose one explicit compatibility strategy:
   - migrate the existing column/data to `password_hash`; or
   - change the model temporarily to match the legacy schema, only if its password storage is secure and compatible.
4. Prepare and verify an external backup plus an exact migration/rollback plan.
5. Obtain separate approval, execute the bounded migration, and rerun the focused login test against a disposable copy of the real schema.
6. Perform a real browser login and verify an authenticated page, logout, and failed-password behavior.
7. Check `SESSION_COOKIE_SECURE = True` for local HTTP. Even after authentication succeeds, a Secure session cookie may not persist over `http://127.0.0.1`; configure this by environment rather than weakening production defaults.

## Handoff cautions

- Do not run `db.create_all()` against the shipped database as a substitute for a migration. It will not safely reconcile incompatible existing columns.
- Do not rename or overwrite `database.db` without the destructive-action approval gate.
- Do not delete the untracked `instance/` directory without previewing its exact contents and obtaining approval.
- Do not claim the README credentials work until a real browser login reaches an authenticated page.
