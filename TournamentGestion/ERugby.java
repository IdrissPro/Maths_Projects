public class ERugby extends EFoot{

    public ERugby(String name) {
        super(name);
    }

    @Override
    int calculNbPoints(Match match){
        int res = super.calculNbPoints(match);
        if (match.getScoreEquipe() >= 30){
            res += 1;
        }
        if (match.getScoreAdversaire() <= 10){
            res += 1;
        }
        return res;
    }
}
