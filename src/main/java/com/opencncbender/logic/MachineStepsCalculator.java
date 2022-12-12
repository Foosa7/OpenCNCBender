package com.opencncbender.logic;

import java.awt.geom.Point2D;


/**
 * Created by Piec on 2019-12-22.
 */

/*
This code is a Java method that calculates a sequence of steps for a machine to perform a specific bending operation on a wire. The method takes an input object of type BendingSteps and returns an object of type Steps that represents the calculated steps for the machine.

The method begins by creating a new Steps object called machineSteps that will be used to store the calculated steps. It then declares several local variables that will be used throughout the method to perform calculations and store intermediate values. These variables include:

angleWithOverbend, which represents the angle of the bend with an added overbend.
arcLength, which represents the length of the arc formed by the bend.
wireDistance, which represents the distance the wire needs to be moved before being bent.
wireRotationAngle, which represents the angle the wire needs to be rotated before being bent.
wireBendAngle, which represents the angle at which the wire will be bent.
lastExtension, which represents the extension of the wire from the previous bending step.
currentExtension, which represents the extension of the wire for the current bending step.
bendDirection, which represents the direction in which the wire will be bent (0 for clockwise, 1 for counterclockwise).
currentStep, which is a SingleStep object that represents the current bending step being processed.
currentMachineStep, which is a SingleMachineStep object that represents the current step calculated for the machine.
The method then enters a loop that iterates over each SingleStep object in the input BendingSteps object. For each step, the method performs several calculations to determine the required distance, rotation, and bend angle for the machine to perform the step. These values are then used to create a new SingleMachineStep object, which is added to the machineSteps object.

Finally, the method returns the machineSteps object that contains the calculated steps for the machine to perform the bending operation.
*/


public class MachineStepsCalculator {

    //private Steps<SingleStep> bendingSteps;
    private MachineGeometry machineGeometry;
    private WireParameters wireParameters;

    private Steps<SingleMachineStep> machineSteps;

    private double minimumArcRadius;

    //private ArrayList<SingleStep> machineStepsList;

    public MachineStepsCalculator(MachineGeometry machineGeometry, WireParameters wireParameters) {
        //this.bendingSteps = bendingSteps;
        this.machineGeometry = machineGeometry;
        this.wireParameters = wireParameters;
        this.minimumArcRadius = machineGeometry.getPinRadius() + wireParameters.getDiameter()/2;

        //machineStepsList = new ArrayList<>(bendingSteps.size());
    }

    public Steps calculateMachineSteps(BendingSteps bendingSteps){

        machineSteps = new Steps();

        double angleWithOverbend, arcLength, wireDistance, wireRotationAngle, wireBendAngle;
        double lastExtension = 0.0;
        double currentExtension = 0.0;
        byte bendDirection;

        SingleStep currentStep;
        SingleMachineStep currentMachineStep;

        //System.out.println(bendingSteps.size());

        for(int i=0; i < bendingSteps.size(); i++){

            currentStep = bendingSteps.get(i);
            angleWithOverbend = calculateAngleWithOverbend(currentStep.getAngleA());
            //System.out.println(angleWithOverbend);

            arcLength = 2 * Math.PI * minimumArcRadius * Math.abs(angleWithOverbend)/360;
            currentExtension = minimumArcRadius * Math.tan(Math.toRadians(Math.abs(angleWithOverbend)/2));

            /*if(currentExtension > currentStep.getDistanceX()/4){
                currentExtension = currentStep.getDistanceX()/4;
            }*/
            //System.out.println(Math.tan(55));
            //System.out.println(minimumArcRadius);
            //System.out.println(arcLength);
            //System.out.println(currentExtension);

            //System.out.println(arcLength);
            //System.out.println(currentExtension);

            //Rotating and moving wire to bend start position
            wireDistance = currentStep.getDistanceX() - lastExtension - currentExtension;
            wireRotationAngle = currentStep.getAngleB();
            currentMachineStep = new SingleMachineStep(wireDistance,wireRotationAngle);
            machineSteps.addStep(currentMachineStep);

            //Bending
            bendDirection = 0;
            if(angleWithOverbend > 0) bendDirection = 1;
            wireBendAngle = calculateMachineAngle(angleWithOverbend);
            currentMachineStep = new SingleMachineStep(arcLength,wireBendAngle,bendDirection);
            machineSteps.addStep(currentMachineStep);

            lastExtension = currentExtension;

            //System.out.println(wireBendAngle);
        }

        return machineSteps;
    }

    private double calculateAngleWithOverbend(double inputAngle){

        double angleWithOverbend;

        if(inputAngle == 0){

            angleWithOverbend = 0;
        }
        else if(inputAngle > 0){

            angleWithOverbend = inputAngle + wireParameters.getOverbendAngle(inputAngle);
        }
        else{

            angleWithOverbend = inputAngle - wireParameters.getOverbendAngle(inputAngle);
        }

        return angleWithOverbend;
    }

    private double calculateMachineAngle(double bendAngle){

        if(bendAngle == 0) {

            return 0;
        }
        else if(bendAngle > 0){

            return calculateMachineAnglePositive(bendAngle);
        }
        else{

            return calculateMachineAngleNegative(bendAngle);
        }
    }

    private double calculateMachineAnglePositive(double bendAngle) {

        Point2D.Double pinCenterPoint, pinTangentPoint, rodCenterPoint;
        Double x,y;

        x = machineGeometry.getPinOffset();
        y = machineGeometry.getPinSpacing()/2;
        pinCenterPoint = new Point2D.Double(x,y);

        x = pinCenterPoint.x + Math.cos(Math.toRadians(bendAngle-90))*(machineGeometry.getPinRadius()+wireParameters.getDiameter());
        y = pinCenterPoint.y + Math.sin(Math.toRadians(bendAngle-90))*(machineGeometry.getPinRadius()+wireParameters.getDiameter());
        pinTangentPoint = new Point2D.Double(x,y);

        if(bendAngle == 90){

            x = pinTangentPoint.x + machineGeometry.getRodRadius();
            y = Math.sqrt(machineGeometry.getBendingRadius()*machineGeometry.getBendingRadius()-x*x);
            rodCenterPoint = new Point2D.Double(x,y);

            return calculateRodAngle(rodCenterPoint);
        }
        else{

            double [] lineParameters; //[m,k]
            Point2D.Double [][] rodPoints;

            lineParameters = calculateLineParameters(bendAngle,pinTangentPoint);
            rodPoints = calculateRodPoints(lineParameters);
            rodCenterPoint = chooseRodCenterPoint(rodPoints,bendAngle);

            //System.out.println(rodCenterPoint.x);
            //System.out.println(rodCenterPoint.y);

            return calculateRodAngle(rodCenterPoint);
        }
    }

    private double calculateMachineAngleNegative(double bendAngle) {

        Point2D.Double pinCenterPoint, pinTangentPoint, rodCenterPoint;
        Double x,y;

        x = machineGeometry.getPinOffset();
        y = -1 * machineGeometry.getPinSpacing()/2;
        pinCenterPoint = new Point2D.Double(x,y);

        x = pinCenterPoint.x + Math.cos(Math.toRadians(90+bendAngle))*(machineGeometry.getPinRadius()+wireParameters.getDiameter());
        y = pinCenterPoint.y + Math.sin(Math.toRadians(90+bendAngle))*(machineGeometry.getPinRadius()+wireParameters.getDiameter());
        pinTangentPoint = new Point2D.Double(x,y);

        if(bendAngle == -90){

            x = pinTangentPoint.x + machineGeometry.getRodRadius();
            y = -1*Math.sqrt(machineGeometry.getBendingRadius()*machineGeometry.getBendingRadius()-x*x);
            rodCenterPoint = new Point2D.Double(x,y);

            return calculateRodAngle(rodCenterPoint);

        }
        else{

            double [] lineParameters; //[m,k]
            Point2D.Double [][] rodPoints;

            lineParameters = calculateLineParameters(bendAngle,pinTangentPoint);
            rodPoints = calculateRodPoints(lineParameters);
            rodCenterPoint = chooseRodCenterPoint(rodPoints,bendAngle);

            return calculateRodAngle(rodCenterPoint);
        }
    }

    private Point2D.Double chooseRodCenterPoint(Point2D.Double[][] rodPoints, double bendAngle) {

        Point2D.Double correctCenterPoint;
        Point2D.Double [] firstSelection = new Point2D.Double[2];
        int k=0;

        if(bendAngle > 0){
            if(bendAngle < 90){

                for(int i=0; i<4; i++){
                    if(rodPoints[i][0].y > rodPoints[i][1].y){
                        firstSelection[k] = rodPoints[i][1];
                        k++;
                    }
                }
                if(firstSelection[0].x > firstSelection[1].x) correctCenterPoint = firstSelection[0];
                else correctCenterPoint = firstSelection[1];
            }
            else{

                for(int i=0; i<4; i++){
                    if(rodPoints[i][0].y < rodPoints[i][1].y){
                        firstSelection[k] = rodPoints[i][1];
                        k++;
                    }
                }
                if(firstSelection[0].x < firstSelection[1].x) correctCenterPoint = firstSelection[0];
                else correctCenterPoint = firstSelection[1];
            }
        }
        else{
            if(bendAngle > -90){

                for(int i=0; i<4; i++){
                    if(rodPoints[i][0].y < rodPoints[i][1].y){
                        firstSelection[k] = rodPoints[i][1];
                        k++;
                    }
                }
                if(firstSelection[0].x > firstSelection[1].x) correctCenterPoint = firstSelection[0];
                else correctCenterPoint = firstSelection[1];
            }
            else{

                for(int i=0; i<4; i++){
                    if(rodPoints[i][0].y > rodPoints[i][1].y){
                        firstSelection[k] = rodPoints[i][1];
                        k++;
                    }
                }
                if(firstSelection[0].x < firstSelection[1].x) correctCenterPoint = firstSelection[0];
                else correctCenterPoint = firstSelection[1];
            }
        }

        return correctCenterPoint;
    }

    private double calculateRodAngle(Point2D.Double rodCenterPoint) {

        /*
        Dot product, signed angle between vector rodCenterPoint and vector [1,0] (X axis), simplified calculations
        */

        return Math.toDegrees(Math.atan2(rodCenterPoint.y,rodCenterPoint.x));
    }

    private double [] calculateLineParameters(double angle, Point2D.Double point){

        double [] lineParameters = new double[2];
        double m,k;

        m = Math.tan(Math.toRadians(angle));
        k = point.y - point.x*m;

        lineParameters[0] = m;
        lineParameters[1] = k;

        return lineParameters;
    }

    private Point2D.Double [][] calculateRodPoints(double [] lineParameters){

        double m,k,s;
        Point2D.Double [][] firstRodPoints;
        Point2D.Double [][] secondRodPoints;

        m = lineParameters[0];
        k = lineParameters[1];

        s = machineGeometry.getRodRadius()*Math.sqrt(m*m+1)+k;
        firstRodPoints = solveEquationsForRodPoints(m,k,s);

        s = k - machineGeometry.getRodRadius()*Math.sqrt(m*m+1);
        secondRodPoints = solveEquationsForRodPoints(m,k,s);

        Point2D.Double [][] rodPoints = new Point2D.Double[4][];
        System.arraycopy(firstRodPoints,0,rodPoints,0,2);
        System.arraycopy(secondRodPoints,0,rodPoints,2,2);

        return rodPoints;
    }

    private Point2D.Double [][] solveEquationsForRodPoints(double m, double k, double s) {

        Point2D.Double [][] rodPoints = new Point2D.Double[2][2];
        double a,b,c,x,y,m1,k1;
        double [] root;
        Point2D.Double rodCenterPoint, rodTangentPoint;

        a = 1 + m*m;
        b = 2*m*s;
        c = s*s - machineGeometry.getBendingRadius()*machineGeometry.getBendingRadius();
        root = quadraticEquation(a,b,c);

        for(int i=0; i < 2; i++){

            x = root[i];
            y = m*x+s;

            rodCenterPoint = new Point2D.Double(x,y);
            rodPoints[i][1]= rodCenterPoint;

            //perpendicular line through rodCenterPoint, finding intersection point
            m1 = -1/m;
            k1 = rodCenterPoint.y - rodCenterPoint.x*m1;
            x = (k - k1)/(m1 - m);
            y = m*x + k;

            rodTangentPoint = new Point2D.Double(x,y);
            rodPoints[i][0]= rodTangentPoint;
        }

        return rodPoints;
    }

    private double [] quadraticEquation(double a, double b, double c){

        double [] root = new double[2];
        double delta;

        delta = b*b - 4*a*c;

        if(delta >= 0){

            root[0] = (-1*b - Math.sqrt(delta))/(2*a);
            root[1] = (-1*b + Math.sqrt(delta))/(2*a);
        }

        return root;
    }
}
