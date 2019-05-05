package gov.nih.ncats.molvec.algo;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import static org.junit.Assert.*;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
/**
 * Created by katzelda on 2/26/19.
 */
@RunWith(Parameterized.class)
@Ignore
public class RegressionTest2 {

        private static List<String> IMPROVED;
        @BeforeClass
        public static void clearOutImprovedList(){
                IMPROVED = new ArrayList<>();
        }

        @AfterClass
        public static void printImproved(){
                if(IMPROVED.isEmpty()){
                    System.out.println("nothing improved");
                     return;
                }
                System.out.println("IMPROVED!!!\n=======");
                IMPROVED.stream().forEach(System.out::println);
        }
    private static void addLARGEST_FRAGMENT_FORMULA_CORRECT(List<Object[]> list, File dir){
        //=================================
        //            LARGEST_FRAGMENT_FORMULA_CORRECT  2
        //--------------------------------------

        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_FORMULA_CORRECT, dir, "cas-30286-75-0"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_FORMULA_CORRECT, dir, "cas-4309-70-0"));
    }

    private static void addATOMS_RIGHT_WRONG_BONDS(List<Object[]> list, File dir){
        //=================================
        //            ATOMS_RIGHT_WRONG_BONDS  14
        //--------------------------------------

        list.add(test(RegressionTestIT.Result.ATOMS_RIGHT_WRONG_BONDS, dir, "cas-141977-79-9"));
        list.add(test(RegressionTestIT.Result.ATOMS_RIGHT_WRONG_BONDS, dir, "cas-16051-77-7"));
        list.add(test(RegressionTestIT.Result.ATOMS_RIGHT_WRONG_BONDS, dir, "cas-16469-74-2"));
        list.add(test(RegressionTestIT.Result.ATOMS_RIGHT_WRONG_BONDS, dir, "cas-17090-79-8"));
        list.add(test(RegressionTestIT.Result.ATOMS_RIGHT_WRONG_BONDS, dir, "cas-52-01-7"));
        list.add(test(RegressionTestIT.Result.ATOMS_RIGHT_WRONG_BONDS, dir, "cas-111490-36-9"));
        list.add(test(RegressionTestIT.Result.ATOMS_RIGHT_WRONG_BONDS, dir, "cas-87-33-2"));
        list.add(test(RegressionTestIT.Result.ATOMS_RIGHT_WRONG_BONDS, dir, "cas-64099-44-1"));
        list.add(test(RegressionTestIT.Result.ATOMS_RIGHT_WRONG_BONDS, dir, "cas-73771-04-7"));
        list.add(test(RegressionTestIT.Result.ATOMS_RIGHT_WRONG_BONDS, dir, "cas-115-44-6"));
        list.add(test(RegressionTestIT.Result.ATOMS_RIGHT_WRONG_BONDS, dir, "cas-135548-15-1"));
        list.add(test(RegressionTestIT.Result.ATOMS_RIGHT_WRONG_BONDS, dir, "cas-124-94-7"));
        list.add(test(RegressionTestIT.Result.ATOMS_RIGHT_WRONG_BONDS, dir, "cas-8069-64-5"));
        list.add(test(RegressionTestIT.Result.ATOMS_RIGHT_WRONG_BONDS, dir, "cas-31721-17-2"));
    }

    private static void addINCORRECT(List<Object[]> list, File dir){
        //=================================
        //            INCORRECT  38
        //--------------------------------------

        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-3093-35-4"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-64521-35-3"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-37517-28-5"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-134457-28-6"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-73232-52-7"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-73747-20-3"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-1767-88-0"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-39791-20-3"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-108050-54-0"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-5611-51-8"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-14149-43-0"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-75734-93-9"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-7232-51-1"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-467-18-5"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-78467-68-2"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-116666-63-8"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-120066-54-8"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-74855-17-7"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-135326-22-6"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-53016-31-2"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-143491-57-0"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-84-97-9"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-35834-26-5"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-5716-20-1"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-54063-44-4"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-9014-89-5"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-292634-27-6"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-39295-97-1"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-53179-09-2"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-108894-41-3"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-163252-36-6"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-141702-36-5"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-114899-77-3"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-114-07-8"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-9005-64-5"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-88931-51-5"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-54143-54-3"));
        list.add(test(RegressionTestIT.Result.INCORRECT, dir, "cas-71990-00-6"));
    }

    private static void addLARGEST_FRAGMENT_ATOM_RIGHT_WRONG_BONDS(List<Object[]> list, File dir){
        //=================================
        //            LARGEST_FRAGMENT_ATOM_RIGHT_WRONG_BONDS  5
        //--------------------------------------

        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_ATOM_RIGHT_WRONG_BONDS, dir, "cas-167465-36-3"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_ATOM_RIGHT_WRONG_BONDS, dir, "cas-5714-76-1"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_ATOM_RIGHT_WRONG_BONDS, dir, "cas-39022-39-4"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_ATOM_RIGHT_WRONG_BONDS, dir, "cas-42879-47-0"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_ATOM_RIGHT_WRONG_BONDS, dir, "cas-3963-95-9"));
    }

    private static void addFORMULA_CORRECT(List<Object[]> list, File dir){
        //=================================
        //            FORMULA_CORRECT  6
        //--------------------------------------

        list.add(test(RegressionTestIT.Result.FORMULA_CORRECT, dir, "cas-52443-21-7"));
        list.add(test(RegressionTestIT.Result.FORMULA_CORRECT, dir, "cas-71116-82-0"));
        list.add(test(RegressionTestIT.Result.FORMULA_CORRECT, dir, "cas-1404-08-6"));
        list.add(test(RegressionTestIT.Result.FORMULA_CORRECT, dir, "cas-105618-02-8"));
        list.add(test(RegressionTestIT.Result.FORMULA_CORRECT, dir, "cas-37115-32-5"));
        list.add(test(RegressionTestIT.Result.FORMULA_CORRECT, dir, "cas-36983-69-4"));
    }

    private static void addLARGEST_FRAGMENT_ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY(List<Object[]> list, File dir){
        //=================================
        //            LARGEST_FRAGMENT_ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY  2
        //--------------------------------------

        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY, dir, "cas-438-41-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY, dir, "cas-518-47-8"));
    }

    private static void addATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY(List<Object[]> list, File dir){
        //=================================
        //            ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY  3
        //--------------------------------------

        list.add(test(RegressionTestIT.Result.ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY, dir, "cas-132236-18-1"));
        list.add(test(RegressionTestIT.Result.ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY, dir, "cas-92812-82-3"));
        list.add(test(RegressionTestIT.Result.ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY, dir, "cas-56038-13-2"));
    }

    private static void addCORRECT_STEREO_INSENSITIVE_INCHI(List<Object[]> list, File dir){
        //=================================
        //            CORRECT_STEREO_INSENSITIVE_INCHI  35
        //--------------------------------------

        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-137-66-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-67337-44-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-22888-70-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-77650-95-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-37529-08-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-13364-32-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-223661-25-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-58944-73-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-89194-77-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-850-52-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-136199-02-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-98769-84-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-75358-37-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-32887-01-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-432-60-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-30418-38-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-642-83-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-39544-74-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-26155-31-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-90243-97-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-17332-61-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-58-18-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-162808-62-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-262352-17-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-17243-38-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-2829-19-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-361-37-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-135928-30-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-576-68-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-170861-63-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-77989-60-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-150683-30-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-221877-54-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-123013-22-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-80225-28-1"));
    }

    private static void addLARGEST_FRAGMENT_CORRECT_STEREO_INSENSITIVE_INCHI(List<Object[]> list, File dir){
        //=================================
        //            LARGEST_FRAGMENT_CORRECT_STEREO_INSENSITIVE_INCHI  6
        //--------------------------------------

        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-54376-91-9"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-34061-34-2"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-107452-79-9"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-21462-39-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-57852-57-0"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_STEREO_INSENSITIVE_INCHI, dir, "cas-149042-61-5"));
    }

    private static void addLARGEST_FRAGMENT_CORRECT_FULL_INCHI(List<Object[]> list, File dir){
        //=================================
        //            LARGEST_FRAGMENT_CORRECT_FULL_INCHI  81
        //--------------------------------------

        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-14176-50-2"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-1050-48-2"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-7327-87-9"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-75859-03-9"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-5588-22-7"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-3546-41-6"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-52618-68-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-70024-40-7"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-2748-88-1"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-82752-99-6"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-81737-62-4"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-6168-86-1"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-82640-04-8"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-125-51-9"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-50978-10-4"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-1239-04-9"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-856-87-1"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-148778-32-9"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-549-18-8"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-53251-94-8"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-5560-62-3"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-15537-76-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-156-51-4"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-6272-74-8"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-2441-88-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-57109-90-7"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-22059-60-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-813-93-4"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-635702-64-6"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-224789-15-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-548-68-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-111523-41-2"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-21817-73-2"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-33286-22-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-61177-45-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-74752-07-1"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-317-52-2"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-532-76-3"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-86641-76-1"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-77495-92-2"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-122111-03-9"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-34866-46-1"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-164579-32-2"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-70476-82-3"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-119271-78-2"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-68576-88-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-550-99-2"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-863127-77-9"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-24358-65-4"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-72318-55-9"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-985-13-7"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-130065-61-1"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-124-90-3"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-76448-47-0"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-11018-89-6"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-2324-94-9"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-7082-21-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-183063-72-1"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-379-79-3"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-56-94-0"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-26786-32-3"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-15826-37-6"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-139781-09-2"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-5634-41-3"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-6576-51-8"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-5579-84-0"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-18833-13-1"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-51552-99-9"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-42408-78-6"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-3151-59-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-225367-66-8"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-279253-83-7"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-4255-24-7"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-3239-45-0"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-102670-59-7"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-22204-29-1"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-65473-14-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-172733-42-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-192564-13-9"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-23736-58-5"));
        list.add(test(RegressionTestIT.Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI, dir, "cas-154-68-7"));
    }

    private static void addERROR(List<Object[]> list, File dir){
        //=================================
        //            ERROR  2
        //--------------------------------------

        list.add(test(RegressionTestIT.Result.ERROR, dir, "cas-39464-87-4"));
        list.add(test(RegressionTestIT.Result.ERROR, dir, "cas-56087-11-7"));
    }

    private static void addCORRECT_FULL_INCHI(List<Object[]> list, File dir){
        //=================================
        //            CORRECT_FULL_INCHI  306
        //--------------------------------------

        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-129-57-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-81-13-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-103844-77-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-60148-52-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-7002-65-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-5626-36-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-104485-01-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-7008-18-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-57935-49-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-97068-30-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-616-91-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-73-24-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-103946-15-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-977-79-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-27076-46-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-38373-83-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-73-09-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-182760-06-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-199396-76-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-24671-26-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-55902-02-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-26605-69-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-97519-39-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-73647-73-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-79243-67-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-22494-27-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-79286-77-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-151287-22-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-31980-29-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-14235-86-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-3579-62-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-15599-27-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-75820-08-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-65886-71-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-446-86-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-2804-00-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-41952-52-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-22161-81-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-56488-58-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-144912-63-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-63245-28-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-110690-43-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-7004-98-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-108001-60-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-139314-01-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-75438-57-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-61263-35-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-87-08-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-65569-29-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-3563-14-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-68307-81-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-18694-40-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-23790-08-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-15866-90-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-111-48-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-112966-96-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-96515-73-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-71617-10-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-15301-48-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-56030-54-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-5728-52-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-51022-76-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-171171-42-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-15992-13-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-1300-94-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-13739-02-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-104145-95-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-5579-95-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-51940-78-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-27164-46-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-99156-66-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-85056-47-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-36144-08-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-32195-33-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-57576-44-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-78967-07-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-172927-65-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-456-59-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-114432-13-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-76-74-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-77-04-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-36590-19-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-103775-14-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-35710-57-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-73-49-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-363-24-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-524-84-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-5741-22-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-378-44-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-138890-62-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-55150-67-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-129-51-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-852-42-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-41826-92-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-127308-82-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-2398-95-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-52042-24-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-170858-33-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-63266-93-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-35273-88-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-64-86-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-531-76-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-70696-66-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-1926-49-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-107-35-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-26807-65-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-69479-26-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-50435-25-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-227318-71-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-77372-61-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-67121-76-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-509-67-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-58186-27-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-132418-36-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-6138-56-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-198481-33-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-170729-80-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-51876-99-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-113-73-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-80349-58-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-16759-59-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-55285-35-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-17692-54-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-5714-73-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-59-33-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-52093-21-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-5964-24-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-25526-93-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-84057-84-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-86487-64-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-66211-92-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-3116-76-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-30271-85-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-118-60-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-34097-16-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-188396-77-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-4317-14-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-179602-65-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-738-70-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-1034-82-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-63388-37-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-28911-01-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-3781-28-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-35700-23-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-55162-26-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-1748-43-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-98-75-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-314-35-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-145414-12-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-32710-91-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-65184-10-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-64241-34-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-147362-57-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-3788-16-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-3820-67-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-61563-18-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-56208-01-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-95-27-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-71-00-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-107489-37-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-59227-89-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-55299-11-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-15687-05-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-4682-36-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-161600-01-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-132-20-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-94-36-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-131796-63-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-89778-27-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-141374-81-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-59-67-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-80844-07-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-101-20-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-129029-23-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-6903-79-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-55149-05-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-76053-16-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-112-72-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-124066-33-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-17692-34-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-132100-55-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-137-26-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-79-01-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-94-26-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-89-25-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-134678-17-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-167305-00-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-525-30-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-57-63-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-10596-23-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-7261-97-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-69014-14-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-112984-60-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-86116-60-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-39567-20-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-26864-56-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-87178-42-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-7455-39-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-2056-56-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-7753-60-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-2363-58-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-13221-27-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-521-74-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-63927-95-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-42461-84-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-84203-09-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-54239-37-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-54063-50-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-33817-09-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-179545-77-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-18296-44-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-187164-19-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-148-56-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-481631-45-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-15534-92-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-120-91-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-53-34-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-127-48-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-25827-76-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-15639-50-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-97-59-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-89943-82-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-502-85-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-149556-49-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-80879-63-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-52157-83-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-50-39-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-90139-06-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-66608-04-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-5580-25-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-77519-25-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-76448-31-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-22494-42-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-144035-83-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-99200-09-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-4444-23-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-80614-21-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-171335-80-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-26513-79-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-594-91-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-94-07-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-69-81-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-56741-95-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-742-20-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-83602-05-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-1641-17-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-149488-17-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-93-88-9"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-25092-41-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-89-83-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-423-55-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-3562-99-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-1642-54-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-59643-91-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-71195-56-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-64-65-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-153-18-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-624-49-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-519-30-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-30924-31-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-518-61-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-21925-88-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-75-19-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-536748-46-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-58662-84-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-74129-03-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-2218-68-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-22365-40-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-195532-12-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-829-74-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-494-03-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-36616-52-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-10331-57-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-78415-72-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-54-92-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-60-46-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-86393-37-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-125372-33-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-74103-07-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-73514-87-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-77360-52-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-431-89-0"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-150812-12-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-6223-35-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-126-93-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-119625-78-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-66969-81-1"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-193153-04-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-18464-39-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-15599-44-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-5581-46-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-209799-67-7"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-36980-34-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-54739-19-4"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-126294-30-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-24815-24-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-58-20-8"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-66984-59-6"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-54340-63-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-510-53-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-89622-90-2"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-148-82-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-66148-78-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-85465-82-3"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-548-00-5"));
        list.add(test(RegressionTestIT.Result.CORRECT_FULL_INCHI, dir, "cas-61-78-9"));
    }

    @Parameterized.Parameters(name="{0}")
    public static List<Object[]> getData(){
        File dir = new File(RegressionTest2.class.getResource("/regressionTest/usan").getFile());

        List<Object[]> list = new ArrayList<>();

        addLARGEST_FRAGMENT_FORMULA_CORRECT(list, dir);
        addATOMS_RIGHT_WRONG_BONDS(list, dir);
        addINCORRECT(list, dir);
        addLARGEST_FRAGMENT_ATOM_RIGHT_WRONG_BONDS(list, dir);
        addFORMULA_CORRECT(list, dir);
        addLARGEST_FRAGMENT_ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY(list, dir);
        addATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY(list, dir);
        addCORRECT_STEREO_INSENSITIVE_INCHI(list, dir);
        addLARGEST_FRAGMENT_CORRECT_STEREO_INSENSITIVE_INCHI(list, dir);
        addLARGEST_FRAGMENT_CORRECT_FULL_INCHI(list, dir);
        addERROR(list, dir);
        addCORRECT_FULL_INCHI(list, dir);
        return list;
    }

    private static Object[] test(RegressionTestIT.Result expected, File dir, String fileNameRoot){

        return new Object[]{fileNameRoot, expected, getImageFor(dir, fileNameRoot), getMoleculeFileFor(dir, fileNameRoot)};
    }
    private static File getMoleculeFileFor(File dir, String nameRoot){
        File smi = new File(dir, nameRoot +".smi");
        if(smi.exists()){
            return smi;
        }
        File sdf = new File(dir, nameRoot +".sdf");
        if(!sdf.exists()){
            throw new IllegalStateException(" not found " + sdf.getAbsolutePath());
        }
        return sdf;
    }
    private static File getImageFor(File dir, String nameRoot){
        File png = new File(dir, nameRoot +".png");
        if(png.exists()){
            return png;
        }
        return new File(dir, nameRoot +".tif");
    }

    private RegressionTestIT.Result expectedResult;
    private File inputImage, expectedSmi;

    public RegressionTest2(String ignored, RegressionTestIT.Result expectedResult, File inputImage, File expectedSmi) {
        this.expectedResult = expectedResult;
        this.inputImage = inputImage;
        this.expectedSmi = expectedSmi;
    }

    @Test
    public void test(){

        RegressionTestIT.Result actual = RegressionTestIT.testMolecule(inputImage, expectedSmi, 300).result;
        if(actual.ordinal() < expectedResult.ordinal()){
                IMPROVED.add( inputImage.getName() + "  " + expectedResult + "  => " + actual);
        }
        assertTrue(expectedResult +" but was " + actual, actual.ordinal() <= expectedResult.ordinal());
    }
}
