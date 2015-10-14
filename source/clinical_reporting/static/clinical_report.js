$(function() {
   $("#report").find('.long_line').load(function(e){
       var self = $(this);
       //var elId = target.attr('id');
       if (self.is(".long_line")) {
           self.css("min_width", self.css("width"));
           //alert('The mouse was over'+ elId);
       }
   });
});