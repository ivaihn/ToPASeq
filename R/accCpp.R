accCpp <-
function (AM, nodeNames, index) 
.Primitive(".Call")(accCpp, AM, nodeNames, 
    index)
    
    downstreamCpp <-
function (AM, nodeNames, index) 
.Primitive(".Call")(downstreamCpp, AM, nodeNames, 
    index)

