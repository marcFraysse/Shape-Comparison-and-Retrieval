import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;

import processing.core.*;
	
public class Main {

	/**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        // setup the top-level enclosing frame
        final JFrame frame = new JFrame("Processing Application");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        // setup the panels
        // panel for visualization
        JPanel panel = new JPanel();
        // panel for buttons
        JPanel buttonPanel = new JPanel();

        // create an initialize the Processing visualization class
        final MeshViewer vView = new MeshViewer();
        vView.init();
        vView.draw();


        // setup the buttons
        JButton bAddElement = new JButton("Add Element");
        // tell the button to use the Processing visualization class
        // as its action listener
        /*bAddElement.addActionListener(vView);*/


        // add the visual elements to the panels
        panel.add(vView);
        buttonPanel.add(bAddElement);
        panel.add(buttonPanel);

        // add the panel in the frame
        frame.add(panel);
        //assign a size for the frame
        //reading the size from the Processing visualization class
        frame.setSize(vView.getSize().width, vView.getSize().height + 100);

        //display the frame
        frame.setVisible(true);
    }

	
}
