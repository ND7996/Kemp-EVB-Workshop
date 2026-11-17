// Check the actual full-enzyme workshop page and the saved-data interaction.
const {chromium}=require('playwright'),assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path');
(async()=>{
 const root=path.join(__dirname,'../kemp_evb');
 const browser=await chromium.launch({executablePath:process.env.CHROMIUM_PATH||'/usr/bin/chromium',headless:true,args:['--no-sandbox']});
 try{
  const page=await browser.newPage({viewport:{width:1280,height:900}});const errors=[];page.on('pageerror',e=>errors.push(e.message));
  await page.setContent(fs.readFileSync(path.join(root,'index.html'),'utf8'));
  assert.match(await page.locator('h1').textContent(),/Run EVB for the Kemp enzyme/);
  assert.equal(await page.locator('header .button').first().getAttribute('href'),'#linux');
  const install=await page.locator('#linux').textContent();
  assert.match(install,/run_kemp_evb\.py/);assert.match(install,/--mode profile --bins 18 --open/);assert.match(install,/\.\/EVBOpenMM/);
  assert.match(await page.locator('#results').textContent(),/LRF\/SCAAS/);
  assert.match(await page.locator('#results').textContent(),/5\.25 ps/);
  assert.equal(await page.locator('#results tbody tr').count(),6);
  assert.match(await page.locator('#profile').textContent(),/9\.57/);assert.match(await page.locator('#profile').textContent(),/13\.2/);
  await page.setViewportSize({width:390,height:844});
  assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth>innerWidth),false,'Mobile horizontal overflow');
  const scan=await browser.newPage();scan.on('pageerror',e=>errors.push(e.message));
  await scan.setContent(fs.readFileSync(path.join(root,'viewers/parameter_scan.html'),'utf8'));
  assert.equal(await scan.locator('#paper-barrier-line').getAttribute('data-value'),'13.2');
  const before=await scan.locator('[data-role=selected-curve]').first().getAttribute('points');
  const values=await scan.locator('#choice option').evaluateAll(nodes=>nodes.map(n=>n.value));
  await scan.locator('#choice').selectOption(values[1]);
  assert.notEqual(await scan.locator('[data-role=selected-curve]').first().getAttribute('points'),before);
  assert.equal(await scan.locator('#paper-barrier-line').getAttribute('data-value'),'13.2');
  const full=await browser.newPage();await full.setContent(fs.readFileSync(path.join(root,'viewers/full_enzyme_recorded_profile.html'),'utf8'));
  assert.equal(await full.locator('pre').count(),0,'Profile must not expose the raw metadata block');
  assert.equal(await full.locator('h2').count(),0,'Profile metadata heading removed');
  assert.match(await full.locator('svg').textContent(),/Grey: published Q6 barrier/);
  assert.deepEqual(errors,[]);
  console.log('Passed: full-enzyme commands, Q6/OpenMM comparison, mobile layout, archived full profile and parameter controls.');
 }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1});
