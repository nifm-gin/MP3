function singleClick(canvas) {
// Applies adaptative zoom when clicking on an image
// If the image is not a reference, it is replaced by the 'large' image
var ctx = canvas.getContext('2d');
var currImg = canvas.nextElementSibling;
console.log(currImg);
switch (currImg.classList.contains("zoomed")){

	case false:	// Image is not zoomed
		currImg.classList.add('zoomed');
		// Distinguish references and maps for scaling 
		if (currImg.closest('tr').rowIndex > 1) { // Image is not a Ref because not on first row (row 0 = names)
			const isRef = 0;
			//var refImg = $('#includedContent').find('tr:eq(1)').find('td:eq('+currImg.closest('td').cellIndex+')')[0].lastElementChild;
            var refImg = $('#includedContent').find('tr:eq(1)').find('td:eq('+currImg.closest('td').cellIndex+')')[0].children[2];
            //console.log($('#includedContent').find('tr:eq(1)').find('td:eq('+currImg.closest('td').cellIndex+')')[0].children[2]);
            //console.log(refImg);
		    S = computeScale(isRef);
		    canvas.width = 2*canvas.width;
		    ctx.clearRect(0, 0, canvas.width, canvas.height);
		    ctx.drawImage(currImg, 0, 0, canvas.width/2, canvas.height);
	    	ctx.drawImage(refImg, currImg.width/2, 0, canvas.width/2, canvas.height);
	    	brightAdjustOnOne([currImg, refImg]);

	    } else { // Image is a Ref
	    	const isRef = 1;
	    	S = computeScale(isRef);
	    	brightAdjustOnOne([currImg]);
	    };
	    $(canvas).css({"transform-origin":"0 0 0", "-webkit-transform": "scale("+S+")", "position": "relative", "z-index": "1", "-webkit-transition": ".1s ease-out", "transition": ".1s ease-out", "opacity": "100%"})
	    break;

	case true:
		currImg.classList.remove('zoomed');
		$(canvas).css({"transform-origin":"0 0 0", "-webkit-transform": "scale(1)", "position": "relative", "z-index": "0", "-webkit-transition": ".0s ease-out", "transition": ".0s ease-out"})
		ctx.clearRect(0, 0, canvas.width, canvas.height);
		
		if (currImg.closest('tr').rowIndex > 1) { // Image is not a Ref because not on first row (row 0 = names)
			canvas.width = canvas.width/2;
	    }
	    ctx.drawImage(currImg, 0,0);
	    brightAdjustOnOne([currImg]);
		break;

	}; //end switch   
};

function drawAll() {
	// Iterates through the images loaded and inserts a canvas element before each of them. The image is then
	// drawn in the canvas, which grants access to pixel values ans allows contrast control
	for (var i = 0; i < document.images.length; i++) {
		// Creates canvas element before each image
		canvas = document.createElement('canvas');
		canvas.setAttribute('width', document.images[i].naturalWidth);
		canvas.setAttribute('height', document.images[i].naturalHeight);
		canvas.onclick = function () { singleClick(this)}
		document.images[i].parentNode.insertBefore(canvas,document.images[i]);
		ctx = canvas.getContext('2d', {alpha: false}); //alpha put to false to disable png transparency
		// Draws the image on the canvas
		ctx.drawImage(document.images[i], 0, 0);
		$('#'+document.images[i].id).data('OriValues', ctx.getImageData(0,0,canvas.width, canvas.height));	
	};
	brightAdjust();
}

function computeScale(isRef) { // WIP 
	// Computes scaling factor to apply when zooming
	var fracWidth = 0.3; // Target fraction of width after scale
	// The aim is to have a constant size after scale, no matter the window zoom
	var S = (fracWidth * $('html').innerWidth()) / 100;
	if (S >=1) {
		switch (isRef){
			case 0:
				return S;
				break;
			case 1:
				return S/2;
				break;
		}
	} else {
		return 1;
	}
};

function showMap(e) {
	// loads the table corresponding to user selection
	var Path = "HTML/";
	var mapToShow = document.getElementById("mapToDisplay");
	$("#includedContent").load(Path+mapToShow.value+".html", function() {
		Srcs = getSources();
		promiseOfAllImages(Srcs)
			.then(function () {
				drawAll();
			})		
		} );

	};

var promiseOfAllImages = function (Srcs) {
	// Wait until ALL images are loaded
	return Promise.all(
	Srcs.map(function (s) {
		// Load each tile, and "resolve" when done
		return new Promise(function (resolve) {
			var img = new Image();
			img.src = s;
			img.onload = function () {
				// Image has loaded... resolve the promise!
				resolve(img);
			};
		});
	})    
	);
};

function getSources() {
	// fetches src attribute from all image elements
	var Imgs = document.images; 
	var imgSrcs = [];
	for (var i =0; i<Imgs.length; i++) {
		imgSrcs.push(Imgs[i].src);
	}
	return imgSrcs
}

function pxVariation(imgDataOri, imgData, value) {
	// Computes pixel values after brightness change
	for (var pixel = 0; pixel < imgData.data.length; pixel++) {
		tmp = imgDataOri.data[pixel] + parseInt(value,10);
		if (tmp > 255) {
			imgData.data[pixel] = 255;
		} else if (tmp < 0) {
			imgData.data[pixel] = 0;
		} else {
			imgData.data[pixel] = tmp;
		}
	}
	return imgData
}

function brightAdjust() { 
	// Runs when brightness slider is clicked and when clicking on 'show', applies on all images
	var slider = document.getElementById("brightRange");
	var imgList = document.images;
	var value = slider.value;
	for (var i = 0; i < imgList.length; i++) {
		var canvas = imgList[i].previousElementSibling;
		var ctx = canvas.getContext('2d');
		if (imgList[i].classList.contains('zoomed') && imgList[i].closest('tr').rowIndex > 1) {
			var Ref = $('#includedContent').find('tr:eq(1)').find('td:eq('+imgList[i].closest('td').cellIndex+')')[0].lastElementChild;
			brightAdjustOnOne([imgList[i], Ref]);

		} else {
			brightAdjustOnOne([imgList[i]]);
		}
	};
}

function brightAdjustOnOne(imgList) {
	// Applies brightness to one image (or a pair if ref needed)
	var slider = document.getElementById("brightRange");
	var value = slider.value;
	var canvas = imgList[0].previousElementSibling;
	var ctx = canvas.getContext('2d');
	switch (imgList.length){
		case 1: // Only one image (dezoom)		
			var imgDataOri = $("#"+imgList[0].id).data('OriValues');
			var imgData = ctx.createImageData(imgDataOri.width, imgDataOri.height);
			imgData = pxVariation(imgDataOri, imgData, value);
			ctx.putImageData(imgData, 0, 0);
			break;

		case 2: // Two images (zoom)
				var imgDataOri = $("#"+imgList[0].id).data('OriValues');
				var imgData = ctx.createImageData(imgDataOri.width, imgDataOri.height);
				imgData = pxVariation(imgDataOri, imgData, value);
				var imgDataOriRef = $("#"+imgList[1].id).data('OriValues');
				var imgDataRef = ctx.createImageData(imgDataOriRef.width, imgDataOriRef.height);
				imgDataRef = pxVariation(imgDataOriRef, imgDataRef, value);
				createImageBitmap(imgData).then(function(response) {ctx.drawImage(response, 0, 0, canvas.width/2, canvas.height)});
				createImageBitmap(imgDataRef).then(function(response) {ctx.drawImage(response, canvas.width/2, 0, canvas.width/2, canvas.height)});
				break;		
	};			
}
