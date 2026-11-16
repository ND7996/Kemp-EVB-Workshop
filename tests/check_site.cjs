// Developer check; saved workshop pages themselves need only a browser.
const {chromium}=require('playwright'),assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path');
(async()=>{
 const root=path.join(__dirname,'../kemp_evb');const browser=await chromium.launch({executablePath:process.env.CHROMIUM_PATH||'/usr/bin/chromium',headless:true,args:['--no-sandbox']});
 try{
 const page=await browser.newPage({viewport:{width:1280,height:900}});const errors=[];page.on('pageerror',e=>errors.push(e.message));await page.setContent(fs.readFileSync(path.join(root,'index.html'),'utf8'));
 assert.match(await page.locator('#matrix-result').textContent(),/E− = -10.000; E\+ = 10.000/);
 await page.locator('#matrix-h').fill('0');assert.match(await page.locator('#matrix-result').textContent(),/Degenerate/);
 await page.locator('#matrix-p').fill('20');assert.match(await page.locator('#matrix-result').textContent(),/R weight = 1.000, P weight = 0.000/);
 await page.locator('#matrix-alpha').fill('-30');assert.match(await page.locator('#matrix-result').textContent(),/R weight = 0.000, P weight = 1.000/);
 await page.setViewportSize({width:390,height:844});assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth>innerWidth),false,'Mobile overflow');
 const slides=await browser.newPage();slides.on('pageerror',e=>errors.push(e.message));await slides.setContent(fs.readFileSync(path.join(root,'presentation.html'),'utf8'));assert.equal(await slides.locator('.slide').count(),16);await slides.keyboard.press('End');assert.match(await slides.locator('.slide.active h2').textContent(),/Grey/);
 const scan=await browser.newPage();scan.on('pageerror',e=>errors.push(e.message));await scan.setContent(fs.readFileSync(path.join(root,'viewers/parameter_scan.html'),'utf8'));assert.match(await scan.locator('#baseline-label').textContent(),/not Q6/);assert.match(await scan.locator('#plot').textContent(),/not the authors’ Q6 profile/);
 assert.deepEqual(errors,[]);console.log('Matrix eigenvalues/weights, mobile layout,16-slide navigation and grey-baseline label passed.');
 }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1});
