public class EVolley extends Equipe {
    public EVolley(String name) {
        super(name);
    }

    @Override
    int calculNbPoints(Match match){
        int res;
        res = match.getScoreEquipe();
        return res;
    }
}
