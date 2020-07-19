package business.signalProcessing;

import business.lie.SE2;
import business.lie.Se2;
import business.lie.SO2;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/* Discrete Extended Kalman Filter on Lie Group */

public class DLGEKalmanFilter {

    private SE2 xa, xp;

    private RealMatrix Pa, Pp, Q, R;

    private double time;

    private SO2 W;

    private RealMatrix I;

    public DLGEKalmanFilter(RealVector X, RealVector nr, double time, double dt) {

        /* first measure */
        double r = X.getEntry(0);
        double alpha = X.getEntry(1);
        double vr = X.getEntry(2);

        /* white Gaussian noise */
        double var_r = nr.getEntry(0);
        double var_alpha = nr.getEntry(1);
        double var_vr = nr.getEntry(2);

        double x = r*Math.cos(alpha);
        double y = r*Math.sin(alpha);
        xa = new SE2(x,y,alpha,vr,0);

        double [][] pax = {
                {r,r/dt,0,0,0},
                {r/dt,2*r/dt,0,0,0},
                {0,0,r,r/dt,0},
                {0,0,r/dt,2*r/dt,0},
                {0,0,0,0,r}};
        Pa = MatrixUtils.createRealMatrix(pax).scalarMultiply(10);

        R = MatrixUtils.createRealDiagonalMatrix(new double [] {var_r,var_alpha,var_vr}).scalarMultiply(dt);

        RealMatrix Qx = MatrixUtils.createRealIdentityMatrix(3);
        RealMatrix B = MatrixUtils.createRealMatrix(new double [][] {{0.5*dt*dt,0,0},{dt,0,0},{0,0.5*dt*dt,0},{0,dt,0},{0,0,1}});
        Q = B.multiply(Qx.multiply(B.transpose()));

        I = MatrixUtils.createRealIdentityMatrix(5);

        this.time = time;
    }

    public SE2 updateModel(SE2 s, RealVector q, double dt) {

        RealVector x = s.v();

        double x1 = x.getEntry(0);
        double y1 = x.getEntry(1);
        double theta1 = x.getEntry(2);
        double v1 = x.getEntry(3);
        double w1 = x.getEntry(4);

        double x2;
        double y2;
        if (0!=w1) {
            x2 = x1 + v1 * (Math.cos(theta1) * Math.sin(w1 * dt) / w1 - Math.sin(theta1) * (1 - Math.cos(w1 * dt)) / w1) + q.getEntry(0);
            y2 = y1 + v1 * (Math.sin(theta1) * Math.sin(w1 * dt) / w1 + Math.cos(theta1) * (1 - Math.cos(w1 * dt)) / w1) + q.getEntry(1);
        } else {
            x2 = x1 + v1 * Math.cos(theta1) * dt + q.getEntry(0);
            y2 = y1 + v1 * Math.sin(theta1) * dt + q.getEntry(1);
        }
        double theta2 = theta1 + w1*dt + q.getEntry(2);
        double v2 = v1 + q.getEntry(3);
        double w2 = w1 + q.getEntry(4);

        return new SE2(x2,y2,theta2,v2,w2);
    }

    /* propagation step */
    public void predict(double dt) {

        xp = xa.multiply(new SE2(xa.omega(dt)));

        RealMatrix C = MatrixUtils.createRealMatrix(new double [][] {
                {0,0,0,dt, 0},
                {0,0,0, 0, 0},
                {0,0,0, 0,dt},
                {0,0,0, 0, 0},
                {0,0,0, 0, 0}});
        RealMatrix F = (new SE2(xp.omega(dt).mapMultiply(-1))).Adg().add( Se2.phi(xp.omega(dt)).multiply(C) );

        Pp = F.multiply(Pa).multiply(F.transpose()).add( Se2.phi(xa.omega(dt)).multiply(Q).multiply(Se2.phi(xa.omega(dt)).transpose()) );
    }

    /* update step */
    public RealVector estimate(RealVector X) {

        SO2 Z = new SO2(X);

        RealMatrix H = buildH(xp,X);
        RealMatrix S = H.multiply(Pp).multiply(H.transpose()).add(R);
        RealMatrix K = Pp.multiply(H.transpose()).multiply(MatrixUtils.inverse(S));

        SO2 xpSo2 = xp.h();

        RealVector m = K.operate((new SO2(xpSo2.inv().multiply(Z.g()))).v());

        SE2 mSE2 = new SE2(m);

        xa = xp.multiply(mSE2);
        Pa = Se2.phi(m).multiply(I.subtract(K.multiply(H))).multiply(Pp).multiply(Se2.phi(m).transpose());

        return xa.h().v();
    }

    private RealMatrix buildH(SE2 xp, RealVector X) {

        double r = X.getEntry(0);
        double alpha = X.getEntry(1);
        double vr = X.getEntry(2);

        RealVector xpHat = xp.v();
        double x = xpHat.getEntry(0);
        double y = xpHat.getEntry(1);
        double theta = xpHat.getEntry(2);
        double v = xpHat.getEntry(3);
        double w = xpHat.getEntry(4);

        double sinTheta = Math.sin(theta);
        double cosTheta = Math.cos(theta);
        double r2 = x*x + y*y;
        double r1 = Math.sqrt(x*x + y*y);

        double H11 = ( x*sinTheta - y*cosTheta)/r2;
        double H12 = ( x*cosTheta + y*sinTheta)/r2;
        double H21 = ( x*cosTheta + y*sinTheta)/r1;
        double H22 = (-x*sinTheta + y*cosTheta)/r1;

        double H31 = -v*Math.sin(alpha-theta)*H11;
        double H32 = -v*Math.sin(alpha-theta)*H12;
        double H33 =  v*Math.sin(alpha-theta);
        double H34 =  Math.cos(alpha-theta);

        RealMatrix H = MatrixUtils.createRealMatrix(new double [][] {
                {H11, H12, 0,   0,   0},
                {H21, H22, 0,   0,   0},
                {H31, H32, H33, H34, 0}});
        return H;
    }

}
