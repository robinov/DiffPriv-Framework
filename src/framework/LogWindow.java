package framework;

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.WindowConstants;

public class LogWindow extends JFrame {
	private static final long serialVersionUID = 1L;
	private JTextArea statusLog = new JTextArea(); 
	
	public LogWindow() {
		super("");
		setSize(450,300);
		add(new JScrollPane(statusLog));
		setVisible(true);
		statusLog.setEditable(false);
		setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
	}
	
	public void showInfo(String data) {
		statusLog.append(data);
		this.validate();
	}
}
