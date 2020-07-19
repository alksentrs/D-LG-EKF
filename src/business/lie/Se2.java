package business.lie;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import util.MyMatrixUtils;

/* se2 Lie Algebra */

public class Se2 {

    private RealMatrix g;
    private RealMatrix adg;
    private RealVector nr;
    private double theta;


    public Se2(double x, double y, double theta, double v, double w) {

        RealMatrix aux = MatrixUtils.createRealMatrix(new double [][] { {Math.sin(theta), 1-Math.cos(theta)},
                                                                        {Math.cos(theta)-1, Math.sin(theta)}});

        double [] p = aux.scalarMultiply(theta/(2*(1-Math.cos(theta)))).operate(new double [] {x,y});

        double px = p[0];
        double py = p[1];

        g = MatrixUtils.createRealMatrix(new double [][] {  { 0,     -theta, px, 0, 0, 0},
                                                            { theta,  0,     py, 0, 0, 0},
                                                            { 0,      0,     0,  0, 0, 0},
                                                            { 0,      0,     0,  0, 0, v},
                                                            { 0,      0,     0,  0, 0, w},
                                                            { 0,      0,     0,  0, 0, 0}});

        adg = MatrixUtils.createRealMatrix(new double [][] {    { 0,     -theta,  py, 0, 0},
                                                                { theta,  0,     -px, 0, 0},
                                                                { 0,      0,      0,  0, 0},
                                                                { 0,      0,      0,  0, 0},
                                                                { 0,      0,      0,  0, 0}});

        this.theta = theta;
        nr = MatrixUtils.createRealVector(new double [] {x,y,theta,v,w});
    }

    public Se2(RealVector rv) {

        nr = rv;
        double x = rv.getEntry(0);
        double y = rv.getEntry(1);
        theta = rv.getEntry(2);
        double v = rv.getEntry(3);
        double w = rv.getEntry(4);

        RealMatrix aux = MatrixUtils.createRealMatrix(new double [][] { {Math.sin(theta), 1-Math.cos(theta)},
                                                                      {Math.cos(theta)-1, Math.sin(theta)}});

        double [] p = aux.scalarMultiply(theta/(2*(1-Math.cos(theta)))).operate(new double [] {x,y});

        double px = p[0];
        double py = p[1];

        g = MatrixUtils.createRealMatrix(new double [][] {  { 0,     -theta, px, 0, 0, 0},
                                                            { theta,  0,     py, 0, 0, 0},
                                                            { 0,      0,     0,  0, 0, 0},
                                                            { 0,      0,     0,  0, 0, v},
                                                            { 0,      0,     0,  0, 0, w},
                                                            { 0,      0,     0,  0, 0, 0}});

        adg = MatrixUtils.createRealMatrix(new double [][] {
                { 0,     -theta,  py, 0, 0},
                { theta,  0,     -px, 0, 0},
                { 0,      0,      0,  0, 0},
                { 0,      0,      0,  0, 0},
                { 0,      0,      0,  0, 0}});
    }

    public Se2(RealMatrix g) {
        this.g = g;
        double px = g.getEntry(0,2);
        double py = g.getEntry(1,2);
        theta = g.getEntry(1,0);
        double v = g.getEntry(3,5);
        double w = g.getEntry(4,5);

        adg = MatrixUtils.createRealMatrix(new double [][] {
                { 0,     -theta,  py, 0, 0},
                { theta,  0,     -px, 0, 0},
                { 0,      0,      0,  0, 0},
                { 0,      0,      0,  0, 0},
                { 0,      0,      0,  0, 0}});

        RealMatrix aux = MatrixUtils.createRealMatrix(new double [][] {
                {Math.sin(theta), 1-Math.cos(theta)},
                {Math.cos(theta)-1, Math.sin(theta)}});

        double [] xy = MatrixUtils.inverse(aux).scalarMultiply((2*(1-Math.cos(theta)))/theta).operate(new double [] {px,py});
        double x = xy[0];
        double y = xy[1];

        nr = MatrixUtils.createRealVector(new double [] {x,y,theta,v,w});
    }

    /* ∨ Operator,  x ∈ g -> x ∈ Rn */
    public RealVector v() {
        return nr;
    }

    /* adjoint representation of Rp */
    public Se2 ad(Se2 o) {
        return new Se2(adjoint(o.g()));
    }

    private RealMatrix adjoint(RealMatrix o) {
        return g.multiply(o).subtract(o.multiply(g));
    }

    public RealMatrix adg() {
        return adg;
    }

    /* Exponential Map x ∈ g -> X ∈ G */
    public SE2 exp() {
        return new SE2(MyMatrixUtils.expm(g));
    }

    public Se2 multiply(Se2 o) {
        return new Se2(g.multiply(o.g));
    }

    public Se2 multiply(RealMatrix o) {
        RealMatrix gg = g.multiply(o);
        return new Se2(gg);
    }

    public RealVector operate(RealVector v) {
        return g.operate(v);
    }

    public RealMatrix g() {
        return g;
    }

    public Se2 transpose() {
        return new Se2(g.transpose());
    }

    public static RealMatrix phi(RealVector a) {
        return jacobian(a.mapMultiply(-1));
    }

    /* left Jacobian */
    public static RealMatrix jacobian(RealVector a) {
        double a1 = a.getEntry(0);
        double a2 = a.getEntry(1);
        double a3 = a.getEntry(2);
        double a4 = a.getEntry(3);
        double a5 = a.getEntry(4);

        RealMatrix adg = MatrixUtils.createRealMatrix(new double [][] {
            {0, -a3,  a2, 0, 0},
            {a3,  0, -a1, 0, 0},
            {0,   0,   0, 0, 0},
            {0,   0,   0, 0, 0},
            {0,   0,   0, 0, 0} });

        RealMatrix I = MatrixUtils.createRealIdentityMatrix(a.getDimension());
        RealMatrix j;
        if (a3!=0) {
            j = I.subtract(adg.scalarMultiply((1 - Math.cos(a3)) / (a3 * a3))).add(adg.multiply(adg).scalarMultiply((a3 - Math.sin(a3)) / (a3 * a3 * a3)));
        } else {
            j = I;
        }
        return j;
    }
}
