$(document).ready(function(){
	$(window).bind('scroll',function(){
		var $topPosition = parseInt($(this).scrollTop());
		if($topPosition>=100){
			if(!$('.navBar').hasClass('fixBar')){
				$('.navBar').addClass('fixBar');
			}
				$('.navContent').css('borderBottom','0 solid #54bb98');
				$('.navDetail').stop().slideUp('fast');

			if(!$('.helpGuide').hasClass('fixGuide')){
				$('.helpGuide').addClass('fixGuide');
			}
		}else{
			if($('.navBar').hasClass('fixBar')){
				$('.navBar').removeClass('fixBar');
			}
			$('.navContent').css('borderBottom','1px solid #54bb98');
			$('.navDetail').stop().slideDown('fast');
			if($('.helpGuide').hasClass('fixGuide')){
				$('.helpGuide').removeClass('fixGuide');
			}
		}
	});



	// 点击显现结果列表
	$('.navContent li').bind('click',function(event){
			event.stopPropagation();
			$('.navDetail').stop().slideDown('fast');
			$('.navContent').css('borderBottom','1px solid #54bb98');
	});
	$('.navDetail').bind('mouseleave',function(){
			if($('.helpGuide').hasClass('fixGuide')){
				$('.navDetail').stop().slideUp('fast');
				$('.navContent').css('borderBottom','0px solid #54bb98');
			}
	});
	//隐藏于显示
	(function(){
		var flag=true;
		$('.guideTog').bind('click',function(){
				$(this).prev().toggle('fast');
				$(this).prev().prev().toggle('fast');
				$(this).children('span').toggleClass('spanC');	
				if($(this).children('span').hasClass('spanC')){
					if(!$('.helpGuide').hasClass('helpGuideH')){
						$('.helpGuide').addClass('helpGuideH');
					}	
				}else{
					$('.helpGuide').removeClass('helpGuideH');
				}
				if(flag){
					$('.contentWrap').css('width','1060px');
				}else{
					$('.contentWrap').css('width','840px');
				}
				flag=(!flag);
		});
		$('.guideTog').trigger('click');
	})();
	// 图片墙
		$('.smallImgList ul li').bind('mouseenter',function(){
			var imgSrc = $(this).children('img').attr('src');
			var $currentImgBox = $(this).parents('.imgList').children('.currentImg').find('img');
			$currentImgBox.attr('src',imgSrc);
			$currentImgBox.parent().attr('href',imgSrc);
		});
		$('.toNext').bind('click',function(){
			var $imgsLi = $(this).parent().children('.smallImgList').find('li');
			var $imgsUl = $(this).parent().children('.smallImgList').find('ul');
			if($imgsLi.length>3){
				if(-(parseInt($imgsUl.css('marginTop'))/100)<$imgsLi.length-3){
					$imgsUl.css('marginTop',parseInt($imgsUl.css('marginTop'))-100);
				}
			}
		});
		$('.toPrev').bind('click',function(){
			var $imgsUl = $(this).parent().children('.smallImgList').find('ul');
			if(parseInt($imgsUl.css('marginTop'))/100<0){
					$imgsUl.css('marginTop',parseInt($imgsUl.css('marginTop'))+100);			
			}
		});
	//导航滚动
	(function(){
		if($('.detailWrap ul').length<5){
			$('.nav-go').hide('fast');
		}


		$('.goRight').bind('click',function(){
			var l = $('.detailWrap ul').length;
			if(l>4){
				var currentP = parseInt($(this).prev().css('marginLeft'));
				if(currentP>(4-l)*180){
					$('.goLeft').show('fast');
					$(this).prev().css('marginLeft',currentP-180);
					if(currentP==(5-l)*180){
						$(this).hide('fast');
					}
				}
			}
		});



		$('.goLeft').bind('click',function(){
			var l = $('.detailWrap ul').length;
			if(l>4){
				var currentP = parseInt($(this).next().css('marginLeft'));
				if(currentP<0){
					$('.goRight').show('fast');
					$(this).next().css('marginLeft',currentP+180);
					if(currentP == -180){
						$(this).hide('fast');
					}
				}
			}
		});
	})();
});