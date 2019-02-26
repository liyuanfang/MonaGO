
var app = angular.module('MonaGO', []);
if (size == 0){
    var go_inf = content["go_inf"]
    var array_order = content["array_order"]
    var clusterHierData = content["clusterHierData"]
    var goNodes = content["goNodes"]
    var matrix = content["matrix"]
}
app.constant('go_inf', go_inf);

app.controller('searchCtrl', ['$scope', 'go_inf', function($scope,go_inf) {
    $scope.goArr = [];

/*    for (var property in go_inf) {
	    if (go_inf.hasOwnProperty(property)) {

	        $scope.go_inf.push({"GO_name":property,"GO_id":goNodes[property]['GO_id']});
	    }
	}	*/
    var uniqueContainer = [];

    go_inf.forEach(function(d,i){
        $scope.goArr.push({"id":d.GO_id,"name":d.GO_name});
        d.genes.forEach(function(d){
            if(!uniqueContainer.includes(d)){
                uniqueContainer.push(d);
                $scope.goArr.push({"gene":d});
            }
        })
        
    });

    //}
    
    $scope.isNotEmpty = function(searchText){

        if (searchText.length!=0){
            return true;
        }
        else{
            return false;
        }
    }

    $scope.setSearchBox = function(GO){
        if(GO.gene){
            $scope.searchText = GO.gene;
        }
        else{
            $scope.searchText = GO.id;
        }
        $('#searchBox').remove();

        setTimeout(function(){
            var e = jQuery.Event("keydown");
            e.which = 13; 
            $('#filter').trigger(e);
        }, 100);
 
    }

}]);
