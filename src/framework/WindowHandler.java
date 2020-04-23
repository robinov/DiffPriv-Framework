package framework;

import java.util.logging.Formatter;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogManager;
import java.util.logging.LogRecord;
import java.util.logging.SimpleFormatter;

public class WindowHandler extends Handler {
	private LogWindow window;
	private Formatter formatter;
	private static WindowHandler handler;
	
	public WindowHandler() {
		formatter = new SimpleFormatter();
		LogManager manager = LogManager.getLogManager();
		String className = this.getClass().getName();
		String level = manager.getProperty(className + ".level");
		setLevel(level != null ? Level.parse(level) : Level.INFO);
		if (window == null) {
			window = new LogWindow();
		}
	}
	
	public static synchronized WindowHandler getInstance() {
		if (handler == null) {
			handler = new WindowHandler();
		}
		return handler;
	}
	
	public synchronized void publish(LogRecord record) {
		String message = null;
		if (!isLoggable(record)) { return; }
		message = formatter.format(record);
		window.showInfo(message);
	}
	
	public void close() {
	}
	
	public void flush() {
	}
}
