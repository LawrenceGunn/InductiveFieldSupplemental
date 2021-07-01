classdef (ConstructOnLoad) RectangularCircuitsCommon
    methods (Static)
        
        function wires = wiresFromDimensions(circuitDimensions, offset)
            wires = RectangularCircuitsCommon.wiresFromCorners(...
                {...
                    [circuitDimensions(1)/2, -circuitDimensions(2)/2,0.0],...
                    [circuitDimensions(1)/2, circuitDimensions(2)/2,0.0],...
                    [-circuitDimensions(1)/2, circuitDimensions(2)/2,0.0],...
                    [-circuitDimensions(1)/2, -circuitDimensions(2)/2, 0.0]...
                },...
                offset,...
                true...
            );
        end
 
        function wireList = wiresFromCorners(wireCorners, offsetCornersBy, isClosed)
            numCorners = length(wireCorners);
            numWires = numCorners-1;
            if isClosed
                numWires = numWires+1;
            end
            wireList = cell(numWires,1);

            function corner = addCorner(startIndex, endIndex)
                corner.start = wireCorners{startIndex} + offsetCornersBy;
                corner.end = wireCorners{endIndex} + offsetCornersBy;
                corner.wireId = startIndex;
            end
            
            for iCorner = 1:numCorners-1
                wireList{iCorner} = addCorner(iCorner, iCorner + 1);
            end
            if isClosed
                wireList{numWires} = addCorner(numCorners, 1);
            end
        end
        
    end
end
                       