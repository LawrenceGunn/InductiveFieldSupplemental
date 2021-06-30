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
            numWires = length(wireCorners);
            if isClosed
                numWires = numWires+1;
            end
            wireList = cell(numWires);

            function corner = addCorner(startIndex, endIndex)
                corner.start = wireCorners(startIndex) + offsetCornersBy;
                corner.end = wireCorners(endIndex) + offsetCornersBy;
                corner.wireId = startIndex;
            end
            
            for iCorner = 1:length(wireCorners)
                wireList(iCorner) = addCorner(iCorner, iCorner + 1);
            end
            if isClosed
                wireList(numWires) = addCorner(iCorner, 1);
            end
        end
        
    end
end
                       