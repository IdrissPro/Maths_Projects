public enum Sport {
    FOOT("Equipe de football"),
    FUTSAL("Equipe de Futsal"),
    HAND("Equipe de Hand-ball"),
    VOLLEY("Equipe de Volley-ball"),
    RUGBY("Equipe de Rugby");

    public String getLabel() {
        return label;
    }

    private final String label;

    Sport(String label) {
        this.label = label;
    }

}

