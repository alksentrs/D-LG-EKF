package business.lie;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import util.MyMatrixUtils;

/* so2 Lie Algebra */

public class So2 {

    private RealMatrix g;

    public So2(RealMatrix g) {
        this.g = g;
    }

    /* Exponential Map x ∈ g' -> X ∈ G' */
    public SO2 exp() {
        return new SO2(MyMatrixUtils.expm(g));
    }


}
