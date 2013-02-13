package nl.esciencecenter.estars.glue;

import javax.swing.JFrame;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class AmuseLib {

    private final AmuseFrame frame;

    public AmuseLib() {

        frame = new AmuseFrame();

        frame.setVisible(true);

        frame.appendText("hurray");
    }

    public void addScene(Scene scene) {
        frame.appendText("adding scene at time " + scene.getTime());
    }
}
