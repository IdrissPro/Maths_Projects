public class EHand extends Equipe{

    public EHand(String name) {
        super(name);
    }

    @Override
    int calculNbPoints(Match match){
        int res;
        if (match.getScoreAdversaire() > match.getScoreEquipe()){
            res = 0;
        }
        else if (match.getScoreEquipe() == match.getScoreAdversaire()){
            res = 1;
        }
        else {
            res = 2;
        }
        return res;
    }
}

