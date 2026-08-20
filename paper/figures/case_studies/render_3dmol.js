const fs = require("fs");
const path = require("path");
const { chromium } = require("playwright");

async function main() {
  if (process.argv.length !== 3) {
    throw new Error("Usage: node render_3dmol.js <config.json>");
  }

  const configPath = path.resolve(process.argv[2]);
  const config = JSON.parse(fs.readFileSync(configPath, "utf8"));
  const scriptPath = path.resolve(config.threeDmolScript);
  const outputPath = path.resolve(config.output);
  fs.mkdirSync(path.dirname(outputPath), { recursive: true });

  const models = config.models.map((model) => ({
    ...model,
    pdbText: fs.readFileSync(path.resolve(model.path), "utf8"),
  }));

  const browserCandidates = [
    config.browserExecutable,
    "C:/Program Files/Google/Chrome/Application/chrome.exe",
    "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe",
    "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe",
    "C:/Program Files/Microsoft/Edge/Application/msedge.exe",
  ].filter(Boolean);
  const browserExecutable = browserCandidates.find((candidate) => fs.existsSync(candidate));
  const browser = await chromium.launch({
    headless: true,
    ...(browserExecutable ? { executablePath: browserExecutable } : {}),
    args: [
      "--use-angle=swiftshader",
      "--enable-webgl",
      "--ignore-gpu-blocklist",
      "--disable-dev-shm-usage",
    ],
  });
  const page = await browser.newPage({
    viewport: { width: config.width || 1500, height: config.height || 1000 },
    deviceScaleFactor: config.deviceScaleFactor || 2,
  });
  await page.setContent(`
    <!doctype html>
    <html>
      <head><meta charset="utf-8"></head>
      <body style="margin:0; overflow:hidden; background:#ffffff;">
        <div id="viewer" style="position:absolute; inset:0;"></div>
      </body>
    </html>
  `);
  await page.addScriptTag({ path: scriptPath });

  await page.evaluate(async ({ config, models }) => {
    const viewer = $3Dmol.createViewer("viewer", {
      backgroundColor: config.background || "white",
      antialias: true,
    });

    const loaded = [];
    for (const item of models) {
      const model = viewer.addModel(item.pdbText, "pdb");
      loaded.push(model);
      const baseSelection = item.selection || {};
      const baseStyle = item.baseStyle || { cartoon: { color: "#D9D9D9" } };
      model.setStyle(baseSelection, baseStyle);

      for (const entry of item.residueStyles || []) {
        const selection = {
          ...(entry.selection || {}),
        };
        model.setStyle(selection, entry.style || {});
      }
      for (const entry of item.addStyles || []) {
        model.setStyle(entry.selection || {}, entry.style || {});
      }
    }

    for (const line of config.lines || []) {
      viewer.addLine({
        start: line.start,
        end: line.end,
        color: line.color || "#555555",
        dashed: line.dashed !== false,
        linewidth: line.linewidth || 2,
      });
    }
    for (const label of config.labels || []) {
      viewer.addLabel(label.text, {
        position: label.position,
        font: "Arial",
        fontSize: label.fontSize || 18,
        fontColor: label.fontColor || "#222222",
        backgroundColor: label.backgroundColor || "#FFFFFF",
        backgroundOpacity: label.backgroundOpacity ?? 0.82,
        borderColor: label.borderColor || "#DDDDDD",
        borderThickness: label.borderThickness ?? 1,
        inFront: true,
      });
    }

    if (config.focusSelection) {
      viewer.zoomTo(config.focusSelection);
    } else {
      viewer.zoomTo();
    }
    if (config.zoom) viewer.zoom(config.zoom);
    for (const rotation of config.rotations || []) {
      viewer.rotate(rotation.angle, rotation.axis);
    }
    viewer.render();
    await new Promise((resolve) => setTimeout(resolve, config.renderDelayMs || 900));
  }, { config, models });

  await page.screenshot({ path: outputPath, omitBackground: false });
  await browser.close();
}

main().catch((error) => {
  process.stderr.write(`${error.stack || error}\n`);
  process.exit(1);
});
