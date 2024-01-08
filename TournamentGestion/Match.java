public class Match {
    private int scoreEquipe;

    private int scoreAdversaire;

    public void majScore(int scoreEquipe, int scoreAdversaire) {
        this.scoreEquipe = scoreEquipe;
        this.scoreAdversaire = scoreAdversaire;
    }

    public int getScoreEquipe() {
        return scoreEquipe;
    }

    public int getScoreAdversaire() {
        return scoreAdversaire;
    }

}
