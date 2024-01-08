public abstract class Equipe {
    private final String name;
    private int nbMatch;
    private int nbButsEncaisses;
    private int nbButsMarques;
    private int pointsEquipe;

    public Equipe(String name) {
        this.name = name;
    }
    abstract int calculNbPoints(Match match);

    private float calculMoyenneEncaisses(){
        return (float) nbButsEncaisses / nbMatch;
    }

    private float calculMoyenneMarques(){
        return (float) nbButsMarques / nbMatch;
    }

    public void majSuiteAUnMatch(Match match){
        pointsEquipe=pointsEquipe + calculNbPoints(match);
        nbMatch=nbMatch+1;
        nbButsEncaisses=nbButsEncaisses+match.getScoreAdversaire();
        nbButsMarques=nbButsMarques+match.getScoreEquipe();
        System.out.println("Match equipe "+name+" score:"+match.getScoreEquipe()+"-"+match.getScoreAdversaire());
    }

    public String restituerResultat(){
        return name+"\t"+"nbPoints : "+pointsEquipe+"\t"+
                "nbButsMarques : "+nbButsMarques+"\t"+
                "moyMarques : "+calculMoyenneMarques()+"\t"+
                "nbButsEncaisses : "+nbButsEncaisses+"\t"+
                "moyEncaissee : "+calculMoyenneEncaisses();
    }

}
