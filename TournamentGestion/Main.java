import java.util.Scanner;
import java.io.*;
public class Main {
    public static void main(String[] args) {

        //On crée d'abord un championnat vide
        Championnat championnat = new Championnat();
        System.out.println("Entrez le nombre d'équipes du club participant au championnat");
        Scanner scan= new Scanner(System.in);
        int nbEquipes=scan.nextInt();

        //on entre les résultats pour nos différentes équipes
        int j;
        for (int i=0;i<nbEquipes;i++)
        {
            j=i+1;
            System.out.println("Sélectionnez le type de sport de l'équipe "+j+": ");
            championnat.entrerScoreUtilisateur(scan);
        }
 
        //l'affichage du championnat est copié dans un fichier
        String resultatCode=championnat.restituerResChampionnat();
        try {
        FileWriter writer = new FileWriter("resultat_championnat.txt");
        writer.write(resultatCode);
        writer.close();
        System.out.println("Le résultat du championnat a aussi été écrit dans un fichier texte.");
        } catch (IOException excep) {
            System.out.println("Une ereur est survenue lors de l'écriture dans le fichier.");
            excep.printStackTrace();
        }

        //enfin, on affiche les résultats de chaque équipe
        System.out.println(resultatCode);
        scan.close();


    }
        
}
