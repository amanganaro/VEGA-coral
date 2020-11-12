/**
 *
 * @author Giovanna Lavado
 * @author Alberto Manganaro
 * 
 */
class StructuralAttribute {
    private String attribute;
    private double corrWeight;

    public StructuralAttribute(String atrribute, double corrWeight) {
        this.attribute = atrribute;
        this.corrWeight = corrWeight;
    }

    public String getAttribute() {
        return attribute;
    }

    public void setAttribute(String attribute) {
        this.attribute = attribute;
    }

    public double getCorrWeight() {
        return corrWeight;
    }

    public void setCorrWeight(double corrWeight) {
        this.corrWeight = corrWeight;
    }

    @Override
    public String toString() {
        return attribute + ": " + corrWeight;
    }
    
}
