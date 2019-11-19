/*-
 * #%L
 * Fiji distribution of ImageJ for the life sciences.
 * %%
 * Copyright (C) 2007 - 2017 Fiji developers.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 * #L%
 */
package utils;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.gui.MultiLineLabel;
import ij.plugin.BrowserLauncher;

import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;


public class GUIHelper
{

    public static final void addHyperLinkListener(final MultiLineLabel text, final String myURL )
    {
        if ( text != null && myURL != null )
        {
            text.addMouseListener( new MouseAdapter()
            {
                @Override
                public void mouseClicked( final MouseEvent e )
                {
                    try
                    {
                        BrowserLauncher.openURL( myURL );
                    }
                    catch ( Exception ex )
                    {
                        IJ.log( "" + ex);
                    }
                }

                @Override
                public void mouseEntered( final MouseEvent e )
                {
                    text.setForeground( Color.BLUE );
                    text.setCursor( new Cursor( Cursor.HAND_CURSOR ) );
                }

                @Override
                public void mouseExited( final MouseEvent e )
                {
                    text.setForeground( Color.BLACK );
                    text.setCursor( new Cursor( Cursor.DEFAULT_CURSOR ) );
                }
            });
        }
    }

}