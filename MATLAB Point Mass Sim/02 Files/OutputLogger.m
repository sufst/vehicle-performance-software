classdef OutputLogger < handle

    properties
        Messages
        BOutput
    end

    methods

        function obj = OutputLogger(bout)
            obj.Messages = cell(3000,1);
            obj.BOutput = bout;
        end

        function addMsg(obj, msg)
            i = find(cellfun('isempty',obj.Messages),1); % Find 1st empty log cell
            if isempty(i)
                obj.Messages = [obj.Messages; cell(3000,1)]; % ran out of space so add some more
                i = find(cellfun('isempty',obj.Messages),1); % Find 1st empty log cell
            end
            obj.Messages{i} = msg;
            if obj.BOutput
                fprintf('%s\n',msg)
            end
        end

        function trim(obj)
            i = find(cellfun('isempty',obj.Messages),1);
            obj.Messages = obj.Messages(1:i-1);
        end
        
    end
end