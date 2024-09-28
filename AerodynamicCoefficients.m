classdef AerodynamicCoefficients
    properties
        materialAspectRatio = 2.25;
        materialLiftSlope = 5.7;
        materialFlapPercentage = 0.2;
        materialZeroLiftAoa = -2.0;
        materialStallAngleHigh = 30.0;
        materialStallAngleLow = -10.0;
        materialSkinFriction = 0.002;
        materialFlapDeflection = 0.0;

        angleOfAttack;
        liftSlope;
        zeroLiftAngleOfAttack;
        stallAngleHigh;
        stallAngleLow;
        liftCoefficient;
        inducedAngle;
        effectiveAngle;
        tangentialCoefficient;
        normalCoefficient;
        dragCoefficient;
        
        curve;
    end
    
    methods
        function obj = AerodynamicCoefficients()
            obj = obj.UpdateCurve();
        end
        
        function obj = UpdateCurve(obj)
            aoaRange = -10:30; 
            obj.curve = zeros(size(aoaRange));
            
            for i = 1:length(aoaRange)
                obj.angleOfAttack = deg2rad(aoaRange(i));
                obj = obj.CalculateCoefficients();
                obj.curve(i) = obj.dragCoefficient;
            end
            
            figure;
            plot(aoaRange, obj.curve);
            xlabel('\color{red} Angle of Attack');
            ylabel('\color{green} CD');
        end
        
        function obj = CalculateCoefficients(obj)
            obj.liftSlope = obj.materialLiftSlope * (obj.materialAspectRatio / ...
                (obj.materialAspectRatio + 2.0 * (obj.materialAspectRatio + 4.0) / (obj.materialAspectRatio + 2.0)));
            
            theta = acos(2.0 * obj.materialFlapPercentage - 1.0);
            flapEffectivenessFactor = 1.0 - (theta - sin(theta)) / pi;
            viscosityFactor = obj.Remap(abs(rad2deg(obj.materialFlapDeflection)), 0.0, 50.0, 0.3, 0.6);
            liftCoefficientIncrement = obj.liftSlope * flapEffectivenessFactor * viscosityFactor * obj.materialFlapDeflection;
            obj.zeroLiftAngleOfAttack = deg2rad(obj.materialZeroLiftAoa) - liftCoefficientIncrement / obj.liftSlope;
            
            maxLiftCoefficientIncrement = liftCoefficientIncrement * obj.Remap(obj.materialFlapPercentage, 0.0, 0.4, 0.3, 0.6);
            maxLiftCoefficientPositive = obj.liftSlope * (deg2rad(obj.materialStallAngleHigh) - deg2rad(obj.materialZeroLiftAoa)) + maxLiftCoefficientIncrement;
            maxLiftCoefficientNegative = obj.liftSlope * (deg2rad(obj.materialStallAngleLow) - deg2rad(obj.materialZeroLiftAoa)) + maxLiftCoefficientIncrement;
            
            obj.stallAngleHigh = obj.zeroLiftAngleOfAttack + maxLiftCoefficientPositive / obj.liftSlope;
            obj.stallAngleLow = obj.zeroLiftAngleOfAttack + maxLiftCoefficientNegative / obj.liftSlope;
            
            if obj.angleOfAttack > obj.stallAngleHigh || obj.angleOfAttack < obj.stallAngleLow
                obj = obj.CalculateHighAngleOfAttackCoefficients();
            else
                obj = obj.CalculateLowAngleOfAttackCoefficients();
            end
        end
        
        function obj = CalculateLowAngleOfAttackCoefficients(obj)
            obj.liftCoefficient = obj.liftSlope * (obj.angleOfAttack - obj.zeroLiftAngleOfAttack);
            obj.inducedAngle = obj.liftCoefficient / (pi * obj.materialAspectRatio);
            obj.effectiveAngle = obj.angleOfAttack - obj.zeroLiftAngleOfAttack - obj.inducedAngle;
            obj.tangentialCoefficient = obj.materialSkinFriction * cos(obj.effectiveAngle);
            obj.normalCoefficient = (obj.liftCoefficient + obj.tangentialCoefficient * sin(obj.effectiveAngle)) / cos(obj.effectiveAngle);
            obj.dragCoefficient = obj.normalCoefficient * sin(obj.effectiveAngle) + obj.tangentialCoefficient * cos(obj.effectiveAngle);
        end
        
        function obj = CalculateHighAngleOfAttackCoefficients(obj)
            friction90degrees = -0.0426 * obj.materialFlapDeflection^2 + 0.021 * obj.materialFlapDeflection + 1.98;
            
            if obj.angleOfAttack > obj.stallAngleHigh
                lowAoaLiftCoefficient = obj.liftSlope * (obj.stallAngleHigh - obj.zeroLiftAngleOfAttack);
            else
                lowAoaLiftCoefficient = obj.liftSlope * (obj.stallAngleLow - obj.zeroLiftAngleOfAttack);
            end
            
            obj.inducedAngle = lowAoaLiftCoefficient / (pi * obj.materialAspectRatio);
            
            if obj.angleOfAttack > obj.stallAngleHigh
                t = interp1([pi/2, obj.stallAngleHigh], [0, 1], obj.angleOfAttack, 'linear', 'extrap');
            else
                t = interp1([-pi/2, obj.stallAngleLow], [0, 1], obj.angleOfAttack, 'linear', 'extrap');
            end
            
            obj.inducedAngle = (1 - t) * obj.inducedAngle;
            obj.effectiveAngle = obj.angleOfAttack - obj.zeroLiftAngleOfAttack - obj.inducedAngle;
            
            obj.normalCoefficient = friction90degrees * sin(obj.effectiveAngle) * ...
                (1.0 / (0.56 + 0.44 * sin(obj.effectiveAngle)) - 0.41 * (1.0 - exp(-17.0 / obj.materialAspectRatio)));
            
            obj.tangentialCoefficient = 0.5 * obj.materialSkinFriction * cos(obj.effectiveAngle);
            obj.liftCoefficient = obj.normalCoefficient * cos(obj.effectiveAngle) - obj.tangentialCoefficient * sin(obj.effectiveAngle);
            obj.dragCoefficient = obj.normalCoefficient * sin(obj.effectiveAngle) + obj.tangentialCoefficient * cos(obj.effectiveAngle);
        end
        
        function result = Remap(~, value, from1, to1, from2, to2)
            result = (value - from1) / (to1 - from1) * (to2 - from2) + from2;
        end
    end
end
