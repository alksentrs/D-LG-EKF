package business.lie;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import util.MyMatrixUtils;

/* SE2 Lie Group */

public class SE2 {

    private RealMatrix G;
    private RealMatrix Adg;
    private RealVector rn;
    private double x, y, theta, v, w;

    public SE2(double x, double y, double theta, double v, double w) {

        G = MatrixUtils.createRealMatrix(new double [][] {
                { Math.cos(theta), -Math.sin(theta), x, 0, 0, 0},
                { Math.sin(theta),  Math.cos(theta), y, 0, 0, 0},
                { 0              ,  0              , 1, 0, 0, 0},
                { 0              ,  0              , 0, 1, 0, v},
                { 0              ,  0              , 0, 0, 1, w},
                { 0              ,  0              , 0, 0, 0, 1}});

        Adg = MatrixUtils.createRealMatrix(new double [][] {
                { Math.cos(theta), -Math.sin(theta),  y, 0, 0},
                { Math.sin(theta),  Math.cos(theta), -x, 0, 0},
                { 0              ,  0              ,  1, 0, 0},
                { 0              ,  0              ,  0, 1, 0},
                { 0              ,  0              ,  0, 0, 1} });

        this.x = x;
        this.y = y;
        this.theta = theta;
        this.v = v;
        this.w = w;

        rn = MatrixUtils.createRealVector(new double [] {x,y,theta,v,w});
    }

    public SE2(RealVector p) {

        x = p.getEntry(0);
        y = p.getEntry(1);
        theta = p.getEntry(2);
        v = p.getEntry(3);
        w = p.getEntry(4);

        G = MatrixUtils.createRealMatrix(new double [][] {
                { Math.cos(theta), -Math.sin(theta), x, 0, 0, 0},
                { Math.sin(theta),  Math.cos(theta), y, 0, 0, 0},
                { 0              ,  0              , 1, 0, 0, 0},
                { 0              ,  0              , 0, 1, 0, v},
                { 0              ,  0              , 0, 0, 1, w},
                { 0              ,  0              , 0, 0, 0, 1}});

        Adg = MatrixUtils.createRealMatrix(new double [][] {
                { Math.cos(theta), -Math.sin(theta),  y, 0, 0},
                { Math.sin(theta),  Math.cos(theta), -x, 0, 0},
                { 0              ,  0              ,  1, 0, 0},
                { 0              ,  0              ,  0, 1, 0},
                { 0              ,  0              ,  0, 0, 1} });

        rn = MatrixUtils.createRealVector(new double [] {x,y,theta,v,w});

    }

    public SE2(RealMatrix G) {
        this.G = G;
        x = G.getEntry(0,2);
        y = G.getEntry(1,2);
        theta = Math.atan2(G.getEntry(1,0),G.getEntry(0,0));
        v = G.getEntry(3,5);
        w = G.getEntry(4,5);

        Adg = MatrixUtils.createRealMatrix(new double [][] {
                { Math.cos(theta), -Math.sin(theta),  y, 0, 0},
                { Math.sin(theta),  Math.cos(theta), -x, 0, 0},
                { 0              ,  0              ,  1, 0, 0},
                { 0              ,  0              ,  0, 1, 0},
                { 0              ,  0              ,  0, 0, 1} });

        rn = MatrixUtils.createRealVector(new double [] {x,y,theta,v,w});
    }

    private RealMatrix Ad(RealMatrix o) {
        return G.multiply(o.multiply(MatrixUtils.inverse(G)));
    }


    /* ∨ Operator,  X ∈ G -> x ∈ Rn */
    public RealVector v() {
        return rn;
    }

    /* Adjoint representation on G in Rp */
    public RealMatrix Adg() {
        return Adg;
    }

    public RealMatrix g() {
        return G;
    }

    public SE2 add(SE2 o) {
        return new SE2(G.add(o.g()));
    }

    /* system model on Lie groups */
    public RealVector omega(double dt) {
        return MatrixUtils.createRealVector(new double [] {v*dt, 0, w*dt, 0, 0});
    }

    public SE2 multiply(SE2 o) {
        return new SE2(g().multiply(o.g()));
    }

    /* Logarithm Map X ∈ G -> x ∈ g */
    public Se2 log() {
        return new Se2(MyMatrixUtils.logm(G));
    }

    /* h: G -> G' */
    public SO2 h() {
        double alpha = Math.atan2(y,x);
        double r = Math.sqrt(x*x + y*y);
        double vr = v * Math.cos(alpha-theta);
        return new SO2(r,alpha,vr);
    }

    public SE2 transpose() {
        return new SE2(G.transpose());
    }

    public SE2 inv() { return new SE2(MatrixUtils.inverse(G)); }
}
