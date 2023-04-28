import { ComponentFixture, TestBed } from '@angular/core/testing';

import { GeneExpressionComponent } from './gene-expression.component';

describe('GeneExpressionComponent', () => {
  let component: GeneExpressionComponent;
  let fixture: ComponentFixture<GeneExpressionComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ GeneExpressionComponent ]
    })
    .compileComponents();

    fixture = TestBed.createComponent(GeneExpressionComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
