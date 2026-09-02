FROM node:24-alpine AS client-build
WORKDIR /build/client
COPY client/package.json client/package-lock.json ./
RUN npm ci
COPY client/ ./
RUN npm run build

FROM python:3.9-slim AS runtime
WORKDIR /app
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
COPY . ./
COPY --from=client-build /build/static/sequence-workspace ./static/sequence-workspace
CMD ["gunicorn", "flask_bio_app:create_app()", "-b", "0.0.0.0:5000"]
