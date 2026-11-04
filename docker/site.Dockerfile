FROM node:22.22.1-bookworm-slim AS dependencies

WORKDIR /opt/pbmc-ageing/site

RUN apt-get update \
    && apt-get install -y --no-install-recommends python3 \
    && rm -rf /var/lib/apt/lists/*

COPY site/package.json site/package-lock.json ./
RUN npm ci --no-audit --no-fund

FROM dependencies AS browser-dependencies
ENV PLAYWRIGHT_BROWSERS_PATH=/ms-playwright
RUN node node_modules/@playwright/test/cli.js install --with-deps chromium \
    && chmod -R a+rX /ms-playwright

FROM dependencies AS runtime
COPY site/ ./
COPY config/studies.json /opt/pbmc-ageing/config/studies.json
COPY scripts/render_report_index.py /opt/pbmc-ageing/scripts/render_report_index.py

# TypeScript and Vite keep small caches under node_modules. Allow the host UID
# used by Make's Docker runs to write those caches without changing source files.
RUN mkdir -p node_modules/.tmp node_modules/.vite \
    && chmod 1777 node_modules/.tmp node_modules/.vite .

CMD ["node", "node_modules/vite/bin/vite.js", "--help"]

FROM browser-dependencies AS browser
COPY site/ ./
COPY config/studies.json /opt/pbmc-ageing/config/studies.json
COPY scripts/render_report_index.py /opt/pbmc-ageing/scripts/render_report_index.py
RUN mkdir -p node_modules/.tmp node_modules/.vite \
    && chmod 1777 node_modules/.tmp node_modules/.vite .
