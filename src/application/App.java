package application;

import business.lie.SE2;
import business.lie.SO2;
import business.signalProcessing.DLGEKalmanFilter;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealVector;

import java.util.Random;

public class App {

    public static void main(String[] args) {

        double dt = 0.01;
        double time = 0;

        double x = 0;
        double y = 0;
        double theta = 0;
        double v = 200;
        double w = 0;

        /* white Gaussian noise */
        RealVector nr = MatrixUtils.createRealVector(new double [] {50,10,0.03}).mapMultiply(dt);

        SE2 xoSe2 = new SE2(x,y,theta,v,w);
        SO2 xaSo2 = xoSe2.h();

        /* filter initialization */
        DLGEKalmanFilter kalmanFilter = new DLGEKalmanFilter(xaSo2.v(),nr,time,dt);

        Random rnd = new Random();

        for (int i=0; i<100; i++) {

            time += dt;

            /* origin */
            xoSe2 = xoSe2.multiply(new SE2(xoSe2.omega(dt)));
            SO2 xoSo2 = xoSe2.h();

            /* measurement */
            double r = xoSo2.v().getEntry(0) + rnd.nextGaussian()*nr.getEntry(0);
            double alpha = xoSo2.v().getEntry(1) + rnd.nextGaussian()*nr.getEntry(1);
            double vr = xoSo2.v().getEntry(2) + rnd.nextGaussian()*nr.getEntry(2);

            /* filter prediction */
            kalmanFilter.predict(dt);

            /* filter estimation */
            RealVector xe = kalmanFilter.estimate(MatrixUtils.createRealVector(new double [] {r,alpha,vr}));

            System.out.println(time+": ");
            //System.out.println(MatrixUtils.createRealVector(new double [] {r*Math.cos(alpha),r*Math.sin(alpha)}));
            System.out.println(MatrixUtils.createRealVector(new double [] {r,alpha,vr}));
            System.out.println(xoSe2.h().v());
            System.out.println(xe);

        }
    }
}
