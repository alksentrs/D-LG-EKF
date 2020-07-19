package business.lie;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import util.MyMatrixUtils;

/* SO2 Lie Group */

public class SO2 {

    private RealMatrix G;
    private RealVector rn;

    public SO2(double r, double alpha, double vr) {

        G = MatrixUtils.createRealMatrix(new double [][] {  { Math.cos(alpha), -Math.sin(alpha), 0, 0, 0},
                                                            { Math.sin(alpha),  Math.cos(alpha), 0, 0, 0},
                                                            { 0              ,  0              , 1, 0, r},
                                                            { 0              ,  0              , 0, 1, vr},
                                                            { 0              ,  0              , 0, 0, 1}});

        rn = MatrixUtils.createRealVector(new double [] {r,alpha,vr});
    }

    public SO2(RealMatrix G) {
        this.G = G;
        double alpha = Math.atan2(G.getEntry(1,0),G.getEntry(0,0));
        double r = G.getEntry(2,4);
        double vr = G.getEntry(3,4);

        rn = MatrixUtils.createRealVector(new double [] {r,alpha,vr});
    }

    public SO2(RealVector X) {

        double r = X.getEntry(0);
        double alpha = X.getEntry(1);
        double vr = X.getEntry(2);

        G = MatrixUtils.createRealMatrix(new double [][] {  { Math.cos(alpha), -Math.sin(alpha), 0, 0, 0},
                { Math.sin(alpha),  Math.cos(alpha), 0, 0, 0},
                { 0              ,  0              , 1, 0, r},
                { 0              ,  0              , 0, 1, vr},
                { 0              ,  0              , 0, 0, 1}});

        rn = MatrixUtils.createRealVector(new double [] {r,alpha,vr});
    }

    public RealMatrix g() {
        return G;
    }

    /* Logarithm Map X ∈ G' -> x ∈ g' */
    public So2 log() {
        return new So2(MyMatrixUtils.logm(G));
    }

    /* ∨ Operator,  X ∈ G' -> x ∈ Rn */
    public RealVector v() {
        return rn;
    }

    public RealMatrix transpose() { return G.transpose(); }

    public RealMatrix inv() { return MatrixUtils.inverse(G); }
}
