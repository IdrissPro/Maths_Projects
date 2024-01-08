public class EFutsal extends EFoot{

    public EFutsal(String name) {
        super(name);
    }

    @Override
    int calculNbPoints(Match match){
        int res = super.calculNbPoints(match);
        if (match.getScoreEquipe() >= 5){
            res += 1;
        }
        return res;
    }
}

